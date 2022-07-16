from PIL import Image
import dearpygui as dpg
import cmot_tof, freespace_green_resonance, narrow_cooling_tweezer_align
from lib.ew import *




outdir=""

SEQ = cmot_tof.make_sequence("seq", 0.0)

def draw(Supersequence s):
    _connections = s.defaults.connections.items()
    T0 = min([seq.min_time() for seq in s.values()])
    T1 = max([seq.max_time() for seq in s.values()])
    color_C = 0
    hk = 3
    linewidth=0.7
    def ycoord(x, conn, k):
        y = x if isinstance(conn, DigitalConnection) else (x + 10) / 20
        return hk * k + 1 + y

        # set up background lines and histories on each connection
        # each history is a list of matplotlib line objects for each sequence
    conn_history = dict()
    for k, (conn_label, conn) in enumerate(reversed(_connections)):
        if isinstance(conn, DigitalConnection):
            P.ax.axhline(
                ycoord(0, conn, k),
                color="0.85", linestyle="-",
                linewidth=linewidth / 2, zorder=100)
            P.ax.axhline(
                ycoord(1, conn, k),
                color="0.85", linestyle="--",
                linewidth=linewidth / 2, zorder=100)
            elif isinstance(conn, AnalogConnection):
                P.ax.axhline(
                    ycoord(-10.0, conn, k),
                    color="0.85", linestyle="--",
                    linewidth=linewidth / 2, zorder=100)
                P.ax.axhline(
                    ycoord(0.0, conn, k),
                    color="0.85", linestyle="-",
                    linewidth=linewidth / 2, zorder=100)
                P.ax.axhline(
                    ycoord(+10.0, conn, k),
                    color="0.85", linestyle="--",
                    linewidth=linewidth / 2, zorder=100)

plot = SEQ.draw_detailed()
print(plot.get_xlim())
print(plot.get_xticklabels())
print(plot.get_xticks())
print(plot.get_ylim())
print(plot.get_yticklabels())
print(plot.get_yticks())
print(plot.get_ylabel())
print(plot.get_xlabel())
