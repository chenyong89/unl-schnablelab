import os
import pandas as pd


def create_df_from_path(path, pattern='*.png', fs=False):
    """
    generate a dataframe for files following specified pattern in a dir

    e.g:
        GenDataFrameFromPath(path, pattern='*.png')
    """
    fnpaths = list(path.glob(pattern))
    df = pd.DataFrame(dict(zip(['fnpath'], [fnpaths])))
    df['dir'] = df['fnpath'].apply(lambda x: x.parent)
    df['fn'] = df['fnpath'].apply(lambda x: x.name)
    if fs:
        df['size'] = df['fnpath'].apply(lambda x: os.path.getsize(x))
    return df
