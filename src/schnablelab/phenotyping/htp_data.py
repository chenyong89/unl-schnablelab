"""
utils to deal with high throughput phenotyping data
"""

import sys
import pandas as pd

from tqdm import tqdm
from PIL import Image
from pathlib import Path
from shutil import copyfile

from schnablelab.apps.tools import create_df_from_path
from schnablelab.apps.base import ActionDispatcher, OptionParser


def main():
    actions = (
        ("list_img_dir", "list image info given dates or samples"),
        ("summary", "Show summary of image info under project root folder"),
        ("extract_images", "extract RGB images from project root folder"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


class ParseProject:
    def __init__(self, prj_dir_name, sm_idx=1, date_idx=2, time_idx=3):
        self.prj_dir_path = Path(prj_dir_name)
        try:
            df = pd.read_csv("%s.idx.csv" % self.prj_dir_path.name)
            df["fnpath"] = df["fnpath"].apply(lambda x: Path(x))
        except FileNotFoundError:
            print("project index file does not exist, creating one...")
            df = create_df_from_path(self.prj_dir_path, pattern="*")
            df = df[df["fnpath"].apply(lambda x: x.is_dir())]
            df["sm"] = df["fn"].apply(lambda x: x.split("_")[sm_idx])
            df["date"] = df["fn"].apply(lambda x: x.split("_")[date_idx])
            df["time"] = df["fn"].apply(lambda x: x.split("_")[time_idx])
            df = df.sort_values(["sm", "date", "time"]).reset_index(drop=True)
            df.to_csv("%s.idx.csv" % self.prj_dir_path.name, index=False)
        finally:
            self.df = df
        self.sm_counts = self.df["sm"].value_counts().sort_index()
        self.date_counts = self.df["date"].value_counts().sort_index()
        self.SMs = self.df["sm"].unique()
        self.dates = self.df["date"].unique()

    def cond_date(self, date_list):
        for date in date_list:
            if date not in self.dates:
                sys.exit("%s not in the date list" % date)
        cond = self.df["date"].isin(date_list)
        return cond

    def cond_sample(self, sample_list):
        for sm in sample_list:
            if sm not in self.SMs:
                sys.exit("%s not in the sample list" % sm)
        cond = self.df["sm"].isin(sample_list)
        return cond

    def yield_info(
        self, samples=None, dates=None, angle=None, backup_angle=None
    ):
        if samples and not dates:
            df = self.df[self.df["sm"].isin(samples)]
        elif not samples and dates:
            df = self.df[self.df["date"].isin(dates)]
        elif samples and dates:
            df = self.df[
                self.df["sm"].isin(samples) & self.df["date"].isin(dates)
            ]
        else:
            df = self.df.copy()
        pbar = tqdm(df.iterrows(), total=df.shape[0])
        for _, row in pbar:
            sm, d, hms = row["sm"], row["date"], row["time"]
            if angle:
                img_fn = row["fnpath"] / ("Vis_SV_%s/0_0_0.png" % angle)
                if not img_fn.exists():
                    if backup_angle:
                        img_fn = row["fnpath"] / (
                            "Vis_SV_%s/0_0_0.png" % backup_angle
                        )
                    else:
                        print(f"{img_fn} does not exist in the project dir!")
                        continue
            else:
                sys.exit("specif at least one vewing angle for RGB images")
            yield sm, d, hms, img_fn


def list_img_dir(args):
    """
    %prog list_img_dir project_root_dir

    list image folders in the project root directory given dates or samples
    """
    p = OptionParser(list_img_dir.__doc__)
    p.add_option(
        "--samples",
        help="specify samples (comma separated without space if"
        " having more than one samples)",
    )
    p.add_option(
        "--dates",
        help="specify dates (comma separated without space if speciy"
        " more than one dates)",
    )
    opts, args = p.parse_args(args)
    if len(args) == 0:
        sys.exit(not p.print_help())
    (project_folder,) = args
    prj = ParseProject(project_folder)
    if opts.samples and not opts.dates:
        samples = opts.samples.split(",")
        cond = prj.cond_sample(samples)
    elif not opts.samples and opts.dates:
        dates = opts.dates.split(",")
        cond = prj.cond_date(dates)
    elif not opts.samples and not opts.dates:
        print("Specify either samples or dates")
    else:
        print("provide either samples or dates")
    print(prj.df[cond][["fn", "sm", "date", "time"]])


def summary(args):
    """
    %prog summary project_root_dir

    Show summary of image info under project_root_dir
    """
    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    (project_folder,) = args

    prj = ParseProject(project_folder)
    print("Summary of samples:")
    for i, j in prj.sm_counts.items():
        print(i, j)
    print("Summary of dates:")
    for i, j in prj.date_counts.items():
        print(i, j)
    print("Angles for RGB images:")
    for angle in prj.df.loc[0, "fnpath"].glob("Vis_*"):
        print(angle.name)


def extract_images(args):
    """
    %prog extract_images project_root_dir

    extract RGB images from project_root_dir
    """
    p = OptionParser(extract_images.__doc__)
    p.add_option(
        "--out_dir", default=".", help="specify the output image directory"
    )
    p.add_option(
        "--samples",
        help="extract particular samples. multiple samples separated"
        " by comma without space",
    )
    p.add_option(
        "--dates",
        help="extract particular dates. multiple dates separated by"
        " comma without space.",
    )
    p.add_option(
        "--angle",
        default="108",
        help="which viewing angle are you going to extract?",
    )
    p.add_option(
        "--backup_angle",
        help="specify an alternative viewing angle for RGB images if"
        " the above angle does not exist.",
    )
    p.add_option(
        "--copy_only",
        default=False,
        action="store_true",
        help="just copy without resizing and converting image format",
    )
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())
    (project_folder,) = args

    out_dir = Path(opts.out_dir)
    if not out_dir.exists():
        out_dir.mkdir()
        print(f"created output directory '{out_dir}'")

    opts.samples = opts.samples.split(",") if opts.samples else opts.samples
    opts.dates = opts.dates.split(",") if opts.dates else opts.dates

    prj = ParseProject(project_folder)
    for sm, d, hms, path_img_fn in prj.yield_info(
        samples=opts.samples,
        dates=opts.dates,
        angle=opts.angle,
        backup_angle=opts.backup_angle,
    ):
        angle_dir_name = path_img_fn.parts[-2]
        dest_fn = "%s_%s_%s_%s.jpg" % (sm, d, hms, angle_dir_name)
        dest = out_dir / dest_fn
        if dest.exists():
            print(f"{dest} already exists, skip!")
        else:
            if opts.copy_only:
                copyfile(path_img_fn, dest)
            else:
                Image.open(path_img_fn).convert("RGB").resize(
                    (1227, 1028)
                ).save(dest)


if __name__ == "__main__":
    main()
