from click import option, group, argument, echo, Path

from oxfel.conversion import longlist_to_ocelot


@group()
def main():
    """Main entrypoint."""
    name = "oxfel"
    echo(name)
    echo("=" * len(name))
    echo("Build OCELOT models from longlists")


@main.command(no_args_is_help=True)
@argument("longlist", nargs=1, type=Path(exists=True, dir_okay=False))
@option("--config")
@argument("outdir", type=Path(exists=True, file_okay=False))
def convert(longlist, config, outdir):
    """Convert longlist to a set of ocelot files"""
    longlist_to_ocelot(longlist, config, outdir)


if __name__ == "__main__":
    main()  # pragma: no cover, pylint: disable=no-value-for-parameter
