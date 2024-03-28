import matplotlib.pyplot as plt
import latdraw
import pandas as pd
from ocelot.cpbd.magnetic_lattice import MagneticLattice

from . import predefined
from .longlist import get_default_component_list, XFELComponentList



def plot_cathode_to_target(target: str, model_type: str = "real") -> tuple[pd.DataFrame, MagneticLattice, plt.Figure]:
    fn = getattr(predefined, f"cat_to_{target}")
    linac = fn(model_type=model_type)

    twiss, mlat = linac.machine_twiss()

    title = f"Cathode to {target.upper()}, {model_type.capitalize()} Optics"

    fig, (_, ax1, ax2, ax3) = latdraw.plot.three_axes_figure(
        linac.get_sequence(), title=title
    )

    ax1.plot(twiss.s, twiss.beta_x, label=r"$\beta_x$")
    ax1.plot(twiss.s, twiss.beta_y, label=r"$\beta_y$")

    ax2.plot(twiss.s, twiss.Dx, label="$D_x$")
    ax2.plot(twiss.s, twiss.Dy, label="$D_y$")

    ax1.legend()
    ax2.legend()

    ax3.plot(twiss.s, twiss.E)

    ax1.set_ylabel(r"$\beta$ / m")
    ax2.set_ylabel("$D$ / m")
    ax3.set_ylabel("$E$ / GeV")

    ax1.set_ylim(0, twiss.beta_x.median() * 10)
    
    return twiss, mlat, fig


def compare_cathode_to_target(target: str, model_type: str = "real") -> tuple[pd.DataFrame, MagneticLattice, plt.Figure]:
    fn = getattr(predefined, f"cat_to_{target}")
    linac = fn(model_type=model_type)

    twiss, mlat = linac.machine_twiss()

    title = f"Cathode to {target.upper()}, {model_type.capitalize()} Optics"

    fig, (_, ax1, ax2, ax3) = latdraw.plot.three_axes_figure(
        linac.get_sequence(), title=title
    )

    ll = get_default_component_list()
    df = ll[f"I1to{target.upper()}"]

    l1, = ax1.plot(twiss.s, twiss.beta_x, label=r"$\beta_x$, OCELOT")
    ax1.plot(df.S, df.BETX, label=r"$\beta_x$, Component List", color=l1.get_color(), linestyle="--")

    l1, = ax1.plot(twiss.s, twiss.beta_y, label=r"$\beta_y$, OCELOT")
    ax1.plot(df.S, df.BETY, label=r"$\beta_y$, Component List", color=l1.get_color(), linestyle="--")

    l1, = ax2.plot(twiss.s, twiss.Dx, label="$D_x$, OCELOT")
    ax2.plot(df.S, df.DX, label=r"$D_x$, Component List", color=l1.get_color(), linestyle="--")
    l1, = ax2.plot(twiss.s, twiss.Dy, label="$D_y$, OCELOT")
    ax2.plot(df.S, df.DY, label=r"$D_y$, Component List", color=l1.get_color(), linestyle="--")

    ax1.legend(ncol=2)
    ax2.legend(ncol=2)

    l1, = ax3.plot(twiss.s, twiss.E, label="OCELOT")
    l1, = ax3.plot(df.S, df.ENERGY, label="Component List", color=l1.get_color(), linestyle="--")

    ax3.legend()

    ax1.set_ylabel(r"$\beta$ / m")
    ax2.set_ylabel("$D$ / m")
    ax3.set_ylabel("$E$ / GeV")

    # A little something to try to stop the huge betas in the dumps
    # ruining the plot.
    ax1.set_ylim(0, twiss.beta_x.median() * 10)

    return twiss, mlat, fig

def compare_ocelot_lattice_component_list(target: str):
    fn = getattr(predefined, f"cat_to_{target}")
    linac = fn(model_type="real")

    clist = get_default_component_list()[f"I1to{target.upper()}"]
    clist_lattice = component_list_to_latdraw_lattice(clist)
    fig, axes = latdraw.subplots_with_lattices(linac.get_sequence(), clist)

