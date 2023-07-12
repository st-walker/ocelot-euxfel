import toml
import pathlib
from click import option, group, argument, echo, Path

from oxfel.conversion import longlist_to_ocelot
import oxfel.accelerator.lattice
from oxfel.longlist import DEFAULT_LONGLIST


@group()
def main():
    """Main entrypoint."""
    name = "oxfel"
    echo(name)
    echo("=" * len(name))
    echo("Build OCELOT models from longlists")


@main.command(no_args_is_help=True)
@option("--config")
@option("--longlist", nargs=1, type=Path(exists=True, dir_okay=False))
@option("--update", is_flag=True)
# @option("--outdir", type=Path(exists=True, file_okay=False))
def update(longlist, config, update):
    """Convert longlist to a set of ocelot files"""
    outdir = pathlib.Path(oxfel.accelerator.lattice.__file__).parent
    longlist_to_ocelot(DEFAULT_LONGLIST, config, outdir)


if __name__ == "__main__":
    main()  # pragma: no cover, pylint: disable=no-value-for-parameter
