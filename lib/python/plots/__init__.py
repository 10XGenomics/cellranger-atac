import numpy as np
import pandas as pd
from scipy import stats, spatial as spatial


def downsample_scatterplot_by_density(points_df, npoints, dim1, dim2):
    """Return downsampled version of points_df with npoints rows.
    Expects points_df to have columns dim1 and dim2. Downsampling
    is done for these two dimensions.
    Downsampling favors sparse regions in the 2D space formed by
    dim1 and dim2.
    """

    if len(points_df) <= npoints:
        return points_df

    # transform to z(log2) scale so that neighborhood radiuses are well
    # behaved. Dynamic range of both axes will be made about (-1, 1).
    dim1_log2 = np.log2(points_df[dim1] + 1)
    dim1_z = stats.zscore(dim1_log2)
    dim2_log2 = np.log2(points_df[dim2] + 1)
    dim2_z = stats.zscore(dim2_log2)

    # round to two decimals -- max 40k points in (-1,1) 2d grid.
    round_df = pd.DataFrame({'z1': np.round(dim1_z, 2),
                             'z2': np.round(dim2_z, 2)},
                            index=points_df.index)

    # sample from the dups to obtain npoints in total.
    np.random.seed(0)
    is_dup = round_df.duplicated()
    ind_unique = round_df.index[is_dup == False]
    ind_dup = round_df.index[is_dup == True]
    if len(ind_unique) <= npoints:
        samp_dups = np.random.choice(ind_dup,
                                      size=npoints - len(ind_unique),
                                      replace=False)
        return pd.concat([points_df.loc[ind_unique], points_df.loc[samp_dups]])

    # random sampling favors sparse regions, penalizes dense regions.
    tree = spatial.KDTree(round_df.loc[ind_unique])
    radius = 0.1
    neighbors = tree.query_ball_tree(tree, radius)
    frequency = np.array(map(len, neighbors))
    inv_density = radius ** 2 / frequency

    samp_index = np.random.choice(round_df.loc[ind_unique].index,
                                  size=npoints,
                                  replace=False,
                                  p=inv_density / sum(inv_density))

    return points_df.loc[samp_index]
