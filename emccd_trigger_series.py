from lib import *
from lib import CONNECTIONS as C
import time

Nframes = 100
t0 = 0.0 # starting time
dt = 0.4 # spacing between trigger pulses
tau = 1e-3 # width of trigger pulse

comp = MAIN.connect()

indef = False
if indef:
    seq = Sequence.digital_pulse(*C.andor_trig, 0.0, tau)
    comp.enqueue(seq)
    while True:
        comp.rerun()
        time.sleep(dt)
else:
    for k in range(Nframes):
        seq = Sequence.digital_pulse(*C.andor_trig, 0.0, tau)
        comp.enqueue(seq).run().clear()
        time.sleep(dt)
    # seq = Sequence.joinall(*[
    #     Sequence.digital_pulse(*C.andor_trig, t0 + k * dt, tau)
    #     for k in range(Nframes)
    # ])
    # comp.enqueue(seq).run()

comp.clear().disconnect()
