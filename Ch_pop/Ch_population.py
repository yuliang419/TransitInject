import numpy as np
import pandas as pd

dir = ''


# input file
# randomly picked the period bin in the log space.
def random_power(a, b, g, size=1):
    """Power-law gen for pdf(x)\propto x^{g-1} for a<=x<=b"""
    r = np.random.random(size=size)
    ag, bg = a ** g, b ** g
    return (ag + (bg - ag) * r) ** (1. / g)


def pick_planet_radius(m_dwarf=True):
    ## pick_planet_radius
    ## input: currently required no input. Updates are required to cooperate the new input parameter including star
    # types.
    ## output:
    #   - P: period of the planet (day)
    #   - radius: radius of the planet (Earth Radii)

    df = pd.read_csv(dir + 'data/dressing+charbonneau_2015_fine_planet_occurrence_rp_per.txt', \
                     names=['Rp_low', 'Rp_high', 'Rp_mid', 'Per_low', 'Per_hi', 'Per_mid', 'Occ'], skiprows=2)

    df2 = pd.read_csv(dir + 'data/fressin_2013_planet_occurrence_rp_per.txt', \
                      names=['Rp_low', 'Rp_high', 'Rp_mid', 'Per_low', 'Per_hi', 'Per_mid', 'Occ'], skiprows=3)

    # Dressing+ 2015
    if m_dwarf:
        dff = df
    else:
        # Fressin+ 2013
        dff = df2

    prob_per = dff.groupby('Per_mid')['Occ'].sum().values / np.sum(dff.groupby('Per_mid')['Occ'].sum().values)
    Per_mid = dff.groupby('Per_mid')['Per_mid'].mean()
    P_mid = np.random.choice(Per_mid, p=prob_per)
    target_P = dff.loc[np.argmin(np.abs(P_mid - dff.Per_mid))]
    P = 10 ** (np.random.uniform(low=np.log10(target_P.Per_low), high=np.log10(target_P.Per_hi)))

    dfff = dff[(P < dff.Per_hi) & (P > dff.Per_low)]
    prob_rad = dfff.Occ / np.sum(dfff.Occ)
    target = dfff.loc[np.random.choice(dfff.index, p=prob_rad.values)]
    if target.Rp_mid < 1.25:
        radius = np.random.uniform(low=target.Rp_low, high=target.Rp_high)
    else:
        radius = random_power(target.Rp_low, target.Rp_high, -0.7, size=1)[0]
    # -----

    return P, radius

# pick_planet_radius()
# (np.random.uniform(low=np.min(dff.Per_mid),high=np.max(dff.Per_mid)))
