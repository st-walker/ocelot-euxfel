import pandas as pd
from importlib_resources import files

DEFAULT_LONGLIST = files("oxfel.accelerator.lattice") / "component_list_2023.07.01.xls"

def make_default_longlist():
    return XFELLongList(DEFAULT_LONGLIST)


class XFELLongList:
    def __init__(self, xls_path, **read_excel_kwargs):
        kwgs = {"sheet_name": "LONGLIST", "skiprows": [1]}
        kwgs.update(read_excel_kwargs)
        df = pd.read_excel(xls_path, **kwgs)
        self.df = df

    def __getitem__(self, name12):
        try:
            return self.df.set_index("NAME1").loc[name12]
        except KeyError:
            return self.df.set_index("NAME2").loc[name12]

    def name2_to_name1(self, name2):
        df = self.df
        return df[df.NAME2 == name2].NAME1

    def name1_to_name2(self, name1):
        ldf = self.df
        return ldf[ldf.NAME1 == name1].NAME2

    def quadrupole_k1s(self) -> dict[str, dict[str:float]]:
        quads = self.df[self.df["CLASS"] == "QUAD"]
        quad_k1ls = quads["STRENGTH"]
        lengths = quads["LENGTH"]
        name2s = quads["NAME1"]

        quad_k1s = quad_k1ls / lengths

        result = {}
        for name2, quad_k1 in zip(name2s, quad_k1s):
            result[name2] = {"k1": quad_k1}

        return result

    def cavity_settings(self):
        cavities = self.df[self.df["CLASS"] == "LCAV"]
        # Filter TDSs, just pick actual accelerating cavities (and also AH1)
        cavities = cavities[(cavities.TYPE == "C") | (cavities.TYPE == "C3")]

        voltages = cavities["STRENGTH"] * 1e-3  # to GV
        # Phases are in units of 2pi (360).
        # Want htem in degrees for OCELOT...
        phases = cavities["E1/LAG"] * 360
        names = cavities["NAME1"]
        result = {}

        for name, voltage, phase in zip(names, voltages, phases):
            result[name] = {"v": voltage, "phi": phase}
        return result

    def chicane_settings(self):
        bc0 = ["BB.96.I1", "BB.98.I1", "BB.100.I1", "BB.101.I1"]
        bc1 = ["BB.182.B1", "BB.191.B1", "BB.193.B1", "BB.202.B1"]
        bc2 = ["BB.393.B2", "BB.402.B2", "BB.404.B2", "BB.413.B2"]

        bc_dipoles = self.df.set_index("NAME1").loc[bc0 + bc1 + bc2]

        # Rename so they can be accessed in the tuple below
        bc_dipoles["E1"] = bc_dipoles["E1/LAG"]
        bc_dipoles["E2"] = bc_dipoles["E2/FREQ"]

        result = {}
        for tup in bc_dipoles.itertuples():
            result[tup.Index] = {"angle": tup.STRENGTH, "e1": tup.E1, "e2": tup.E2}
        return result

    def focussing_settings(self):
        return self.cavity_settings() | self.quadrupole_k1s() | self.chicane_settings()

    def get_optics_constraint(self, name1) -> dict[str : dict[str, float]]:
        row = self[name1]
        return {
            name1: {
                "beta_x": row.BETX,
                "beta_y": row.BETY,
                "alpha_x": row.ALFX,
                "alpha_y": row.ALFY,
            }
        }
