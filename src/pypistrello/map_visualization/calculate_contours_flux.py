from scipy.ndimage import gaussian_filter
import matplotlib.patheffects as pe

import numpy as np

# Copies the map, fills NaN values with zeros, and applies Gaussian smoothing.
def prepare_map_for_contours(
        image,
        smooth_sigma=1.0,
        fill_nan=False):
    """
    Prepare a 2D image for contour extraction.

    Parameters
    ----------
    image : ndarray
    smooth_sigma : float
    fill_nan : bool

    Returns
    -------
    processed_image
    valid_mask
    """

    mask=np.isfinite(image)

    filled=image.copy()
    filled[~mask]=0

    weight=mask.astype(float)

    smooth_image=gaussian_filter(filled, smooth_sigma)
    smooth_weight=gaussian_filter(weight, smooth_sigma)

    image=smooth_image/smooth_weight

    image[~mask]=np.nan

    return image, mask

def robust_sigma(data):
    data = data[np.isfinite(data)]
    med = np.median(data)
    mad = np.median(np.abs(data-med))
    return 1.4826*mad

# Computes contour levels based on the specified mode and parameters.
def compute_contour_levels(image, params):

    finite = image[np.isfinite(image)]
    mode = params.get("contours_mode", "sigma")

    print(f"Contour mode: {mode}")

    if mode=="manual":
        return np.asarray(params.get("manual_levels"))
    
    elif mode=="percentile":
        return np.percentile(finite,params.get("percentile_levels"))

    elif mode=="sigma":
        sigma = robust_sigma(finite)
        sigma_levels = params.get("sigma_levels", [3, 5, 10, 20])

        return sigma * np.asarray(sigma_levels)

    else:
        raise ValueError("Unknown contour mode")

    
# Save the contours
def save_processed_map(filename, image, levels):

    np.savez_compressed(
        filename,
        image=image,
        levels=levels,
        origin="lower",
        shape=image.shape,
    )

# draw contours on the provided axis

def draw_contours(
        ax,
        image,
        levels,
        contours_outlined=False,
        colors="black",
        linewidths=1.0,
        outline_color="white",
        outline_factor=2.2,
        **kwargs):
    """
    Draw contours, optionally with a white outline.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    image : ndarray
    levels : array-like
    contours_outlined : bool, optional
        If True, draw a thicker white contour underneath.
    colors : str
        Main contour colour.
    linewidths : float
        Width of the main contour.
    outline_color : str
        Colour of the outline.
    outline_factor : float
        Factor by which the outline is thicker than the main contour.
    """

    if contours_outlined:
        ax.contour(
            image,
            levels=levels,
            colors=outline_color,
            linewidths=linewidths * outline_factor,
            **kwargs
        )

    return ax.contour(
        image,
        levels=levels,
        colors=colors,
        linewidths=linewidths,
        **kwargs
    )

# load the contours from a saved file
def load_processed_map(filename):
    data=np.load(filename)

    return (
        data["image"],
        data["levels"]
    )

