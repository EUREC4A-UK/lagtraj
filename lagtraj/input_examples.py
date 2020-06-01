from pathlib import Path
from collections import OrderedDict as OD
import yaml

from asciitree import LeftAligned
from asciitree import drawing


P_ROOT = Path(__file__).parent.parent / "input_examples"


class LagtrajExampleDoesNotExist(Exception):
    pass


INPUT_TYPE_PLURAL = dict(
    trajectory="trajectories", domain="domains", forcing="forcings"
)


def attempt_load(input_name, input_type):
    input_type = INPUT_TYPE_PLURAL.get(input_type, input_type)

    file_path = P_ROOT / input_type / (input_name + ".yaml")
    if not file_path.exists():
        raise LagtrajExampleDoesNotExist

    with open(file_path) as fh:
        defn = yaml.load(fh, Loader=yaml.FullLoader)
        return defn


def print_available(input_types="all"):
    def _print_node(p):
        dirs = p.glob("*")
        if p == P_ROOT and input_types != "all":
            dirs = filter(lambda p_: p_.name in input_types, dirs)

        return OD(
            [(p_node.name.replace(".yaml", ""), _print_node(p_node)) for p_node in dirs]
        )

    tree = {"lagtraj://": _print_node(P_ROOT)}

    box_tr = LeftAligned(draw=drawing.BoxStyle(gfx=drawing.BOX_LIGHT))
    print(box_tr(tree))


def main():
    print(
        "The following domain/forcing/trajectory definitions are"
        " currently included with lagtraj:"
    )
    print("")

    print_available()

    print("\n")
    print(
        "To use for example the `eurec4a_north_atlantic` domain definition\n"
        "for downloading domain data run lagtraj.domain.download as"
        " follows:\n"
        "\n"
        "    $> python -m lagtraj.domain.download"
        " lagtraj://eurec4a_20191209_12_eul"
    )


if __name__ == "__main__":
    main()
