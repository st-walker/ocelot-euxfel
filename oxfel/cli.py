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
from oxfel.longlist import get_default_component_list
from oxfel.astra import (
    load_reference_0320_10k_distribution,
    load_reference_0320_100k_distribution,
)
from oxfel import cat_to_i1d, cat_to_b1d, cat_to_b2d, cat_to_tld, all_models
from oxfel.optics import print_match_point_analysis, START_SIM
from oxfel.plot import plot_cathode_to_target, compare_cathode_to_target, compare_ocelot_lattice_component_list


@group()
def main():
    """Main entrypoint."""
    name = "oxfel"
    echo(name)
    echo("=" * len(name))
    echo("Build OCELOT models from longlists")


@main.command()
@option("--config")
# @argument("longlist", nargs=1, type=Path(exists=True, dir_okay=False))
def update(config):
    """Convert longlist to a set of ocelot files"""
    outdir = pathlib.Path(oxfel.accelerator.lattice.__file__).parent
    longlist_to_ocelot(get_default_component_list(), config, outdir, match=match)


@main.command()
def match():
    match_real_injector()


@main.command()
@argument("target", nargs=1)
@option("--real", is_flag=True)
@option("--compare", is_flag=True)
def plot(target, real, compare):

    model_type = "design"
    if real:
        model_type = "real"

    if compare:
        twiss, _, _ = compare_cathode_to_target(target, model_type=model_type)
    else:
        twiss, _, _ = plot_cathode_to_target(target, model_type=model_type)

    print_match_point_analysis(twiss)

    plt.show()


@main.command()
@argument("target", nargs=1)
def lattice(target):
    fig = compare_ocelot_lattice_component_list(target)
    


@main.command()
def sizes():
    parray032 = load_reference_0320_10k_distribution()
    for model_name, model in all_models(model_type="tracking").items():
        destination = model_name.split("cat_to_")
        print(f"Generating beam size data for {model_name}")
        update_bunch_size_data(destination, model, parray032)



if __name__ == "__main__":
    main()  # pragma: no cover, pylint: disable=no-value-for-parameter
