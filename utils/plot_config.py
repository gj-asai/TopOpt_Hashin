import matplotlib as mpl
import matplotlib.pyplot as plt
import niceplots


def config() -> None:
    plt.style.use(niceplots.get_style())

    mpl.rcParams['axes.spines.right'] = True
    mpl.rcParams['axes.spines.top'] = True

    mpl.rcParams['xtick.direction'] = 'in'
    mpl.rcParams['ytick.direction'] = 'in'

    mpl.rcParams['xtick.top'] = True
    mpl.rcParams['ytick.right'] = True