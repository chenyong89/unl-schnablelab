import cv2
import scipy
import skimage
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.patches import Circle
from PIL import Image


class RGB:
    """
    collection of functions for processing RGB images
    """

    def __init__(self, filename):
        self.fn = filename
        self.format = filename.split(".")[-1]
        self.PIL_img = Image.open(filename)
        self.width, self.height = self.PIL_img.size
        self.array = np.array(self.PIL_img)[:, :, 0:3]
        self.array_r = self.array[:, :, 0]
        self.array_g = self.array[:, :, 1]
        self.array_b = self.array[:, :, 2]
        self.green_idx = (2 * self.array_g) / (
            self.array_r + self.array_b + 0.01
        )
        self.img_dim = self.array.shape

    def crop(self, crp_dim):
        """
        crop_dim: (left, upper, right, lower)
        """
        return self.PIL_img.crop(crp_dim)

    def resize(self, resize_dim):
        """
        resize_dim: (width, height)
        """
        return self.PIL_img.resize(resize_dim)

    def biseg_by_green(self, cutoff=130):
        """
        perform binary segmentation by values in green channel
            convert background pixels to 0 (black)
            convert plant pixels to 255 (white)
        """
        _, thresh = cv2.threshold(self.array_g, cutoff, 255, cv2.THRESH_BINARY)
        thresh_ivt = skimage.util.invert(thresh).astype("uint8")
        return thresh_ivt

    def biseg_by_green_idx(self, cutoff=1.12):
        """
        perform binary segmentation by green index
            convert background pixels to 255 (white)
            convert plant pixels to 0 (black)
        """
        _, thresh = cv2.threshold(
            self.array_gidx, cutoff, 255, cv2.THRESH_BINARY
        )
        thresh = thresh.astype("uint8")
        return thresh

    def get_plant_boundry(self, box=None):
        """
        get boundary of plant from a binary segmentated image
        args:
            box: tuple (left, upper, right, bottom)
        """
        thresh_ivt = self.biseg_by_green()
        if box is not None:
            left, upper, right, bottom = box
            thresh_ivt[0:upper] = 0
            thresh_ivt[bottom:] = 0
            thresh_ivt[:, 0:left] = 0
            thresh_ivt[:, right:] = 0

        height = np.sum(thresh_ivt, axis=1)
        width = np.sum(thresh_ivt, axis=0)
        idx_top = next(x for x, val in enumerate(height) if val > 0)
        idx_bottom = self.height - next(
            x for x, val in enumerate(height[::-1]) if val > 0
        )
        idx_left = next(x for x, val in enumerate(width) if val > 0)
        idx_right = self.width - next(
            x for x, val in enumerate(width[::-1]) if val > 0
        )
        return (idx_top, idx_bottom, idx_left, idx_right)

    def get_convex_hull(self):
        """
        compute the convex hull image of plant
        value of plant pixels: 2
        value of hull region excluding plant: 1 (gray)
        value of non-hull region: 0 (black)
        """
        thresh_ivt = self.biseg_by_green()
        chull = skimage.morphology.convex_hull_image(thresh_ivt)
        chull_diff = np.where(thresh_ivt == 255, 2, chull)
        return chull_diff

    def png_to_jpg(self):
        """
        convert image in png format to jpg format
        """
        return self.PIL_img.convert("RGB")


def call_height(img_fn, zoom_level, threshold=1.12):
    """
    calculate plant height using green index from an RGB image
    """
    if zoom_level == 1:
        upper, bottom, left, right = 60, 1750, 500, 2250
        ratio = 149 / 1925
    if zoom_level == 2:
        upper, bottom, left, right = 180, 1700, 850, 1770
        ratio = 149 / 965
    img = cv2.imread(img_fn)
    img = cv2.cvtColor(img, cv2.COLOR_BGR2RGB)
    img_float = img.astype(np.float)
    img_green_idx = (2 * img_float[:, :, 1]) / (
        img_float[:, :, 0] + img_float[:, :, 2]
    )
    thresh1 = np.where(img_green_idx > threshold, img_green_idx, 0)
    print("remove the chamber border")
    thresh1[0:upper] = 0
    thresh1[bottom:] = 0
    thresh1[:, 0:left] = 0
    thresh1[:, right:] = 0
    thresh1 = (thresh1 / float(thresh1.max())) * 255
    blur = cv2.GaussianBlur(thresh1, (7, 7), 0)
    blur_int = blur.astype(np.uint8)
    _, thresh2 = cv2.threshold(blur_int, 1, 255, cv2.THRESH_BINARY)
    _, contours, _ = cv2.findContours(
        thresh2, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE
    )
    cv2.drawContours(img, contours, -1, (0, 255, 0), 3)

    min_y, max_y = [], []
    for i in contours:
        min_y.append(np.min(i[:, :, 1]))
        max_y.append(np.max(i[:, :, 1]))
    if min_y and max_y:
        y_lowest, y_highest = min(min_y), max(max_y)
        height_pixels = y_highest - y_lowest
        height_cm = height_pixels * ratio
        cv2.line(img, (500, y_lowest), (2000, y_lowest), (255, 0, 0), 7)
        img = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
        new_fn = img_fn.replace(".png", ".height.png")
        print("image with height marked has been saved to {new_fn}")
        cv2.imwrite(new_fn, img)
    return height_cm


def call_organ(rgb_arr, part="stem"):
    crp_shape2d = rgb_arr.shape[0:2]
    if part == "stem":
        r, g, b = 251, 129, 14
    elif part == "panicle":
        r, g, b = 126, 94, 169
    elif part == "leaf":
        r, g, b = 0, 147, 0
    else:
        raise Exception("only support stem, panicle, and leaf")
    p1 = np.full(crp_shape2d, r)
    p2 = np.full(crp_shape2d, g)
    p3 = np.full(crp_shape2d, b)
    p123 = np.stack([p1, p2, p3], axis=2)
    pRGB = np.where(rgb_arr == p123, rgb_arr, 255)
    return pRGB


def color_organ(arr2d, part="stem"):
    cond_k = arr2d == 0
    if part == "stem":
        r, g, b = 251, 129, 14
    elif part == "panicle":
        r, g, b = 126, 94, 169
    elif part == "leaf":
        r, g, b = 0, 147, 0
    else:
        raise Exception("only support stem, panicle, and leaf")
    pr = np.where(cond_k, r, 255)
    pg = np.where(cond_k, g, 255)
    pb = np.where(cond_k, b, 255)
    pRGB = np.stack([pr, pg, pb], axis=2)
    return pRGB


def filter_pixel(arr3d, d=0):
    rgb_img = Image.fromarray(arr3d)
    gray_img = rgb_img.convert(mode="L")
    gray_blur_arr = cv2.GaussianBlur(np.array(gray_img), (3, 3), 0)
    cutoff = (
        pd.Series(gray_blur_arr.flatten())
        .value_counts()
        .index.sort_values()[d]
    )
    arr2d = np.where(gray_blur_arr <= cutoff, 0, 255)
    return arr2d


def polish_organ(img_fn, crop_dim=None, blur_degree=4):
    """
    remove noises from each plant organ using opencv blur
    args:
        crop_dim: tuple (left, upper, right, bottom)
    """
    img = Image.open(img_fn)
    if crop_dim:
        img = np.array(img.crop(crop_dim))
    else:
        img = np.array(img)
    stemRGBraw = call_organ(img, "stem")
    stem = filter_pixel(stemRGBraw)
    stemRGB = color_organ(stem, "stem")
    panicleRGBraw = call_organ(img, "panicle")
    panicle = filter_pixel(panicleRGBraw, d=blur_degree)
    panicleRGB = color_organ(panicle, "panicle")
    leafRGBraw = call_organ(img, "leaf")
    leaf = filter_pixel(leafRGBraw, d=blur_degree)
    leafRGB = color_organ(leaf, "leaf")
    spRGB = np.where(stemRGB == 255, panicleRGB, stemRGB)
    splRGB = np.where(spRGB == 255, leafRGB, spRGB)
    output_img_fn = img_fn.replace(".png", ".polished.png")
    scipy.misc.sm.imsave(output_img_fn, splRGB)


def add_marker(img_fn, marker_coordinates, colors):
    """
    add markers to an image
    args:
        img_fn: image filename
        marker_coordinates:
            2D numpy array with shape (n, 2) where n is number of markers
        colors:
            list of colors for each marker
    return:
        matplotlib ax with markers added
    """
    img = skimage.io.imread(img_fn)
    _, ax = plt.subplots(1)
    ax.imshow(img)
    for xy, c in zip(marker_coordinates, colors):
        circle = Circle(xy, radium=2, color=c)
        ax.add_patch(circle)
    return ax
