import numpy as np
import pandas as pd

path = "observation_data/"

dfDES = pd.read_csv(path + "DESY1.csv")
dfDES = dfDES.assign(
    ylow=dfDES["Y"] - dfDES["-DeltaY"], yhigh=dfDES["+DeltaY"] - dfDES["Y"]
)
dfDES = dfDES.assign(
    xlow=dfDES["X"] - dfDES["-DeltaX"], xhigh=dfDES["+DeltaX"] - dfDES["X"]
)
yerrDES = np.array([dfDES["ylow"], dfDES["yhigh"]])
xerrDES = np.array([dfDES["xlow"], dfDES["xhigh"]])

yerrDESdimless = np.array(
    [
        dfDES["ylow"] * (dfDES["X"] ** 3) / (2 * (np.pi**2)),
        dfDES["yhigh"] * (dfDES["X"] ** 3) / (2 * (np.pi**2)),
    ]
)

# columns = ["k", "P(k)", "delta P(k)+", "delta P(k)-"]
# dfBOSS = pd.read_csv(path + "BOSS.csv", header=None)
# dfBOSS.columns = columns

# dfBOSS = dfBOSS.assign(
#     ylow=dfBOSS["P(k)"] - dfBOSS["delta P(k)-"],
#     yhigh=dfBOSS["delta P(k)+"] - dfBOSS["P(k)"],
# )
# yerrBOSS = np.array([dfBOSS["ylow"], dfBOSS["yhigh"]])

# yerrBOSSdimless = np.array(
#     [
#         dfBOSS["ylow"] * (dfBOSS["k"] ** 3) / (2 * (np.pi**2)),
#         dfBOSS["yhigh"] * (dfBOSS["k"] ** 3) / (2 * (np.pi**2)),
#     ]
# )

dfBOSS = pd.read_csv(path + "BOSS.csv", skiprows=1)

dfBOSS = dfBOSS.assign(
    ylow=dfBOSS["Y"] - dfBOSS["Y.3"], yhigh=dfBOSS["Y.4"] - dfBOSS["Y"]
)
dfBOSS = dfBOSS.assign(
    xlow=dfBOSS["X"] - dfBOSS["X.1"], xhigh=dfBOSS["X.2"] - dfBOSS["X"]
)
yerrBOSS = np.array([dfBOSS["ylow"], dfBOSS["yhigh"]])
xerrBOSS = np.array([dfBOSS["xlow"], dfBOSS["xhigh"]])

yerrBOSSdimless = np.array(
    [
        dfBOSS["ylow"] * (dfBOSS["X"] ** 3) / (2 * (np.pi**2)),
        dfBOSS["yhigh"] * (dfBOSS["X"] ** 3) / (2 * (np.pi**2)),
    ]
)


dfHERA = pd.read_csv(path + "HeraProjected.csv")
dfEDGES = pd.read_csv(path + "EDGESProjected.csv")

dfHERA = dfHERA.assign(
    ylow=dfHERA["Y"] - dfHERA["-DeltaY"], yhigh=dfHERA["+DeltaY"] - dfHERA["Y"]
)
dfHERA = dfHERA.assign(
    xlow=dfHERA["X"] - dfHERA["-DeltaX"], xhigh=dfHERA["+DeltaX"] - dfHERA["X"]
)

HERA_y_dimfull = dfHERA["Y"] / (dfHERA["X"] ** 3) * (2 * (np.pi**2))

yerrHERA = np.array([dfHERA["ylow"], dfHERA["yhigh"]])
xerrHERA = np.array([dfHERA["xlow"], dfHERA["xhigh"]])

yerrHERAdimfull = np.array(
    [
        dfHERA["ylow"] / (dfHERA["X"] ** 3) * (2 * (np.pi**2)),
        dfHERA["yhigh"] / (dfHERA["X"] ** 3) * (2 * (np.pi**2)),
    ]
)

dfEDGES = dfEDGES.assign(
    ylow=dfEDGES["Y"] - dfEDGES["-DeltaY"], yhigh=dfEDGES["+DeltaY"] - dfEDGES["Y"]
)
dfEDGES = dfEDGES.assign(
    xlow=dfEDGES["X"] - dfEDGES["-DeltaX"], xhigh=dfEDGES["+DeltaX"] - dfEDGES["X"]
)

EDGES_y_dimfull = dfEDGES["Y"] / (dfEDGES["X"] ** 3) * (2 * (np.pi**2))

yerrEDGES = np.array([dfEDGES["ylow"], dfEDGES["yhigh"]])
xerrEDGES = np.array([dfEDGES["xlow"], dfEDGES["xhigh"]])

yerrEDGESdimfull = np.array(
    [
        dfEDGES["ylow"] / (dfEDGES["X"] ** 3) * (2 * (np.pi**2)),
        dfEDGES["yhigh"] / (dfEDGES["X"] ** 3) * (2 * (np.pi**2)),
    ]
)

# Adjust projections to agree with lcdm line
dflcdm = pd.read_csv(
    "output/lambdacdm.dat",
    header=None,
    names=["k", "P(k)"],
    skiprows=4,
    delimiter="\s+",
)

dfEDGES["Y"] = [
    1.05 * dflcdm["P(k)"][121] * (dflcdm["k"][121] ** 3) / (2 * np.pi**2),
    dflcdm["P(k)"][129] * (dflcdm["k"][129] ** 3) / (2 * np.pi**2),
]
dfHERA["Y"] = [
    1.05 * dflcdm["P(k)"][121] * (dflcdm["k"][121] ** 3) / (2 * np.pi**2),
    dflcdm["P(k)"][128] * (dflcdm["k"][128] ** 3) / (2 * np.pi**2),
    dflcdm["P(k)"][130] * (dflcdm["k"][130] ** 3) / (2 * np.pi**2),
]


####PLANCK####
dfPlanckTT = pd.read_csv(path + "PlanckTT.csv", skiprows=1)
dfPlanckEE = pd.read_csv(path + "PlanckEE.csv", skiprows=1)
dfPlanckPP = pd.read_csv(path + "PlanckPP.csv", skiprows=1)

dfPlanckTT = dfPlanckTT.assign(
    ylow=dfPlanckTT["Y"] - dfPlanckTT["Y.3"], yhigh=dfPlanckTT["Y.4"] - dfPlanckTT["Y"]
)
dfPlanckTT = dfPlanckTT.assign(
    xlow=dfPlanckTT["X"] - dfPlanckTT["X.1"], xhigh=dfPlanckTT["X.2"] - dfPlanckTT["X"]
)
yerrTT = np.array([dfPlanckTT["ylow"], dfPlanckTT["yhigh"]])
xerrTT = np.array([dfPlanckTT["xlow"], dfPlanckTT["xhigh"]])

yerrTT_dimless = np.array(
    [
        dfPlanckTT["ylow"] * (dfPlanckTT["X"] ** 3) / (2 * (np.pi**2)),
        dfPlanckTT["yhigh"] * (dfPlanckTT["X"] ** 3) / (2 * (np.pi**2)),
    ]
)

dfPlanckEE = dfPlanckEE.assign(
    ylow=dfPlanckEE["Y"] - dfPlanckEE["Y.3"], yhigh=dfPlanckEE["Y.4"] - dfPlanckEE["Y"]
)
dfPlanckEE = dfPlanckEE.assign(
    xlow=dfPlanckEE["X"] - dfPlanckEE["X.1"], xhigh=dfPlanckEE["X.2"] - dfPlanckEE["X"]
)
yerrEE = np.array([dfPlanckEE["ylow"], dfPlanckEE["yhigh"]])
xerrEE = np.array([dfPlanckEE["xlow"], dfPlanckEE["xhigh"]])

yerrEE_dimless = np.array(
    [
        dfPlanckEE["ylow"] * (dfPlanckEE["X"] ** 3) / (2 * (np.pi**2)),
        dfPlanckEE["yhigh"] * (dfPlanckEE["X"] ** 3) / (2 * (np.pi**2)),
    ]
)

dfPlanckPP = dfPlanckPP.assign(
    ylow=dfPlanckPP["Y"] - dfPlanckPP["Y.3"], yhigh=dfPlanckPP["Y.4"] - dfPlanckPP["Y"]
)
dfPlanckPP = dfPlanckPP.assign(
    xlow=dfPlanckPP["X"] - dfPlanckPP["X.1"], xhigh=dfPlanckPP["X.2"] - dfPlanckPP["X"]
)
yerrPP = np.array([dfPlanckPP["ylow"], dfPlanckPP["yhigh"]])
xerrPP = np.array([dfPlanckPP["xlow"], dfPlanckPP["xhigh"]])

yerrPP_dimless = np.array(
    [
        dfPlanckPP["ylow"] * (dfPlanckPP["X"] ** 3) / (2 * (np.pi**2)),
        dfPlanckPP["yhigh"] * (dfPlanckPP["X"] ** 3) / (2 * (np.pi**2)),
    ]
)
