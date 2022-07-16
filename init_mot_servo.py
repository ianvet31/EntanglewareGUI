from lib import *
from lib import CONNECTIONS as C

seq = (
    Sequence.serial_bits_c(
        C.mot3_coils_sig,
        0.1, # choose t = 0.1 just to be safe; s
        AD5791_INIT, 20, # init code is 20 bits
        AD5791_CTRL, 4, # DAC setting code is 4 bits
        C.mot3_coils_clk, C.mot3_coils_sync,
    )
    + Sequence.serial_bits_c(
        C.mot3_coils_sig,
        0.2, # choose t = 0.2 just to be safe; s
        int(441815), 20, # default setting for blue MOT
        AD5791_DAC, 4, # DAC setting code is 4 bits
        C.mot3_coils_clk, C.mot3_coils_sync,
    )
)

comp = MAIN.connect().clear()
comp.enqueue(seq).run()
comp.clear().disconnect()
