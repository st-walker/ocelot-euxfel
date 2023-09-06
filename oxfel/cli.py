import matplotlib.pyplot as plt
import toml
import pathlib
from click import option, group, argument, echo, Path

from oxfel.conversion import (
    longlist_to_ocelot,
    match_real_injector,
    update_bunch_size_data,
)
import oxfel.accelerator.lattice
from oxfel.longlist import DEFAULT_LONGLIST
from oxfel.astra import (
    load_reference_0320_10k_distribution,
    load_reference_0320_100k_distribution,
)
from oxfel import cat_to_i1d, cat_to_b1d, cat_to_b2d, cat_to_tld, all_models
from oxfel.optics import print_match_point_analysis, START_SIM


@group()
def main():
    """Main entrypoint."""
    name = "oxfel"
    echo(name)
    echo("=" * len(name))
    echo("Build OCELOT models from longlists")


@main.command()
@option("--config")
@option("--longlist", nargs=1, type=Path(exists=True, dir_okay=False))
def update(longlist, config):
    """Convert longlist to a set of ocelot files"""
    outdir = pathlib.Path(oxfel.accelerator.lattice.__file__).parent
    longlist_to_ocelot(DEFAULT_LONGLIST, config, outdir, match=match)


@main.command()
def match():
    match_real_injector()


@main.command()
@argument("target", nargs=1)
@option("--design", is_flag=True)
@option("--real", is_flag=True)
@option("--tracking", is_flag=True)
def plot(target, design, real, tracking):
    import latdraw

    if design:
        model_type = "design"
    elif real:
        model_type = "real"
    elif tracking:
        model_type = "tracking"

    if target == "i1d":
        title = f"Cathode to I1D, {model_type.capitalize()} Optics"
        fel = cat_to_i1d(model_type=model_type)
    elif target == "b1d":
        title = f"Cathode to B1D, {model_type.capitalize()} Optics"
        fel = cat_to_b1d(model_type=model_type)
    elif target == "b2d":
        title = f"Cathode to B2D, {model_type.capitalize()} Optics"
        fel = cat_to_b2d(model_type=model_type)
    elif target == "tld":
        title = f"Cathode to TLD, {model_type.capitalize()} Optics"
        fel = cat_to_tld(model_type=model_type)

    if model_type == "tracking":
        parray032 = load_reference_0320_100k_distribution()
        parrayend, twiss = fel.track_optics(parray032, start=START_SIM, physics=True)
    else:
        twiss, mlat = fel.machine_twiss()

    print_match_point_analysis(twiss)

    _, (_, ax1, ax2, ax3) = latdraw.plot.three_axes_figure(
        fel.get_sequence(), title=title
    )

    ax1.plot(twiss.s, twiss.beta_x, label=r"$\beta_x$")
    ax1.plot(twiss.s, twiss.beta_y, label=r"$\beta_y$")

    ax2.plot(twiss.s, twiss.Dx, label="$D_x$")
    ax2.plot(twiss.s, twiss.Dy, label="$D_y$")

    ax1.legend()
    ax2.legend()

    ax3.plot(twiss.s, twiss.E)

    ax1.set_ylabel("$\\beta$ / m")
    ax2.set_ylabel("$D$ / m")
    ax3.set_ylabel("$E$ / GeV")

    plt.show()

@main.command()
def sizes():
    parray032 = load_reference_0320_10k_distribution()
    for model_name, model in all_models(model_type="tracking").items():
        destination = model_name.split("cat_to_")
        print(f"Generating beam size data for {model_name}")
        update_bunch_size_data(destination, model, parray032)

if __name__ == "__main__":
    main()  # pragma: no cover, pylint: disable=no-value-for-parameter
