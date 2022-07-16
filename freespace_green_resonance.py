from lib import *
from lib import CONNECTIONS as C
import numpy as np
from itertools import product
import sys
import time
import pathlib
import timeit

timestamp = get_timestamp()
print(timestamp)
outdir = DATADIRS.freespace_green.joinpath(timestamp)

comments = """
TOF = 0 ms
"""[1:-1]

# CAMERA OPTIONS
camera_config = {
    "exposure_time": 300,
    "gain": 47.99,
    "gamma": 1.25,
    "black_level": 5.0,
    "bin_size": 1,
}

# GLOBAL PARAMETERS
# use `t_<...>` for definite times
# use `tau_<...>` for durations and relative timings
# use `B_<...>` for MOT coil servo settings
# use `shims_<direction>` for shim coil settings
# use `det_<...>` for detunings
# use `f_<...>` for absolute frequencies
# use `delta_<name>` for small steps in corresponding parameters
# use ALL_CAPS for arrays
# add a comment for any literal numbers used

## MAIN SEQUENCE PARAMETERS
# general
reps = 15 # repeat shots for statistics
clk_freq = 10e6 # serial_bits clock frequency; Hz
take_background = True

# timings (not to scale)
# start | t0 | tau_tof | end
#                      |
#                  <tau_flir>
#                <tau_andor> (tau_img)
t0 = 450e-3 # release blue MOT; s
tau_flux_block = -15e-3 # time relative to t0 to stop atom flux; s
tau_tof = 0e-3 # release time before imaging; s
tau_img = 0.5e-3 # use the probe beams to image the atoms after release; s
tau_flir = 0.0e-3 # Flir camera trigger time relative to end of TOF; s
tau_andor = -0.25e-3 # EMCCD camera trigger time relative to end of TOF; s

# coil settings
B_blue = int(441815) # coil setting for blue MOT
bits_Bset = int(20) # number of bits in a servo setting
bits_DACset = int(4) # number of bits in a DAC-mode setting

# frequency parameters
f0 = 93.8 # CMOT AOM frequency; MHz
DET_PROBE = np.linspace(-1500e-3, 1000e-3, 26) # probe beam detuning for imaging; MHz
det_block = 10 # shift the RF on the MOT beams to decouple the fiber; MHz

# power parameters
p0 = 28.0 # power in both probe beam AOMs; dBm

def make_sequence(name: str) -> SuperSequence:
    SEQ = SuperSequence(
        outdir.joinpath("sequences"),
        name,
        {
            "Flir": (Sequence()
                + Sequence.digital_pulse(
                    *C.flir_trig,
                    t0 + tau_tof + tau_flir,
                    camera_config["exposure_time"] * 1e-6 # convert back us -> s
                )
            ).with_color("C2"),

            "EMCCD": (Sequence()
                + Sequence.digital_pulse(
                    *C.andor_trig,
                    t0 + tau_tof + tau_andor - 27.02e-3, # account for shutter
                    tau_img # exposure time = tau_img
                )
            ).with_color("C2"),

            "Load blue MOT": (Sequence()
                + Sequence.digital_hilo(
                    *C.mot3_blue_aom,
                    0.0, # beams on to load atoms from beginning
                    t0
                )
                # + Sequence.digital_hilo(
                #     *C.mot3_blue_sh,
                #     0.0, # beams on to load atoms from beginning
                #     t0 - 5e-3 # close shutter ahead of AOM
                # )
                + Sequence.digital_pulse(
                    *C.push_aom,
                    t0 + tau_flux_block,
                    tau_tof + tau_img
                )
                # + Sequence.digital_pulse(
                #     *C.push_sh,
                #     t0 + tau_flux_block - 2e-3, # account for shutter time
                #     tau_tof + tau_img
                # )
                + Sequence.digital_pulse(
                    *C.mot2_blue_aom, # remove atoms from push beam path to be sure
                    t0 + tau_flux_block,
                    tau_tof + tau_img
                )
                + Sequence.serial_bits_c(
                    C.mot3_coils_sig,
                    0.0, # make sure the coils are on for loading
                    B_blue, bits_Bset,
                    AD5791_DAC, bits_DACset,
                    C.mot3_coils_clk, C.mot3_coils_sync,
                    clk_freq
                )
                + Sequence.serial_bits_c(
                    C.mot3_coils_sig,
                    t0,
                    0, bits_Bset, # turn off for free-space measurement
                    AD5791_DAC, bits_DACset,
                    C.mot3_coils_clk, C.mot3_coils_sync,
                    clk_freq
                )
                + Sequence.serial_bits_c(
                    C.mot3_coils_sig,
                    t0 + tau_tof + tau_img + 1e-3, # make sure the coils are on for loading
                    B_blue, bits_Bset,
                    AD5791_DAC, bits_DACset,
                    C.mot3_coils_clk, C.mot3_coils_sync,
                    clk_freq
                )
            ).with_color("C0"),

            "Shim coils": (Sequence()
                + [ # set shims to default blue MOT values
                    Event.analog(**c, s=c.default) @ (k * 1e-16)
                    for k, c in enumerate([
                        C.shim_coils_fb,
                        C.shim_coils_lr,
                        C.shim_coils_ud,
                    ])
                ]
                + [
                    Event.analog(**c, s=s) @ (t0 + k * 1e-16)
                    for k, (c, s) in enumerate([
                        (C.shim_coils_fb, 0.0),
                        (C.shim_coils_lr, 0.0),
                        (C.shim_coils_ud, 0.0),
                    ])
                ]
            ).with_color("C8"),

            "Green imaging": (Sequence()
                + Sequence.digital_pulse(
                    *C.raman_green_aom,
                    t0 + tau_tof,
                    tau_img
                )
                # + Sequence.digital_pulse(
                #     *C.raman_green_sh,
                #     t0 + tau_tof - 2e-3, # account for 2 ms shutter delay
                #     tau_img
                # )
            ).with_color("C6"),

            "Scope": (
                Sequence.digital_hilo(
                    *C.scope_trig,
                    t0,
                    t0 + tau_tof
                )
            ).with_color("C7"),

        },
        CONNECTIONS
    )
    SEQ["Sequence"] = ( # dummy sequence; ensures that EW backend stays in-sequence
        Sequence.digital_hilo(
            *C.dummy,
            0.0,
            SEQ.max_time() + 1e-3
        )
    ).with_color("k")

    return SEQ

# background sequence -- timings are hard-coded because they're not important
seq_bkgd = SuperSequence(
    outdir.joinpath("sequences"),
    "background",
    {
        "Sequence": (
            Sequence.digital_hilo(*C.dummy, 0.0, 700e-3)
        ).with_color("k"),

        "Push": (
            Sequence.digital_lohi(*C.push_sh, 0.0, 500e-3)
        ).with_color("C3"),

        "MOT beams": (
            Sequence.digital_lohi(*C.mot3_blue_sh, 0.0, 500e-3)
        ).with_color("C0"),

        "MOT coils": (Sequence()
            + Sequence.digital_hilo(*C.mot3_coils_igbt, 0.0, 500e-3)
            + Sequence.digital_hilo(*C.mot3_coils_onoff, 0.0, 500e-3)
        ).with_color("C1"),

        "Green beams": (Sequence()
            + Sequence.digital_hilo(*C.mot3_green_aom, 0.0, 10e-3)
            # + Sequence.digital_hilo(*C.mot3_green_sh, 0.0, 10e-3)
        ).with_color("C6"),

        "Camera": (Sequence
            .digital_pulse(*C.flir_trig, 70e-3, camera_config["exposure_time"] * 1e-6)
        ).with_color("C2"),

    }, CONNECTIONS)

# SCRIPT CONTROLLER
class NarrowCoolingTweezerAlignment(Controller):
    def precmd(self, *args):
        self.comp = MAIN.connect()
        self.comp.set_defaults().def_digital(*C.raman_green_sh, 1)

        self.cam = FLIR.connect()
        (self.cam
            .configure_capture(**camera_config)
            .configure_trigger()
        )

        self.mog = MOGRF.connect()

        # pre-construct all sequences and mogtables so we don't have to do it
        #   in the main loop
        self.names = list()
        self.params = list()
        self.sequences = list()
        self.ssequences = list()
        for d_probe, rep in product(DET_PROBE, range(reps)):
            name = f"tau-tof={tau_tof:.5e}_probe-det={d_probe:.5f}_{rep}"
            self.names.append(name)
            self.params.append((d_probe, rep))
            sseq = make_sequence(name)
            self.sequences.append(sseq.to_sequence())
            self.ssequences.append(sseq)

    def run_sequence(self, *args):
        for (d_probe, rep), seq in zip(self.params, self.sequences):
            if rep == 0:
                (self.mog
                    .set_frequency(3, f0 + d_probe).set_power(3, p0)
                    .set_frequency(4, f0 + d_probe).set_power(4, p0)
                )
            (self.comp
                .enqueue(seq)
                .run()
                .clear()
            )
        if take_background:
            (self.comp
                .enqueue(seq_bkgd.to_sequence())
                .run()
                .clear()
            )

        self.comp.clear().set_defaults().def_digital(*C.raman_green_sh, 1).disconnect()

        (self.mog
            .set_frequency(1, 90.0).set_power(1, 29.04)
            .set_frequency(3, f0)
            .set_frequency(4, f0)
            .disconnect()
        )

    def run_camera(self, *args):
        self.frames = self.cam.acquire_frames(
            num_frames=len(self.params) + int(take_background),
            timeout=5, # s
            roi=[913, 782, 150, 150]  #[978, 737, 20, 20] # [970, 750, 30, 30]
        )
        # might be important to disconnect in the action rather than postcmd
        #   due to possible bad sharing of resources between threads
        self.cam.disconnect()
        self.cam = None

    def on_error(self, ERR: Exception, *args):
        # try to exit gracefully on error, even if it means this part is
        #   rather dirty
        try:
            self.comp.clear().disconnect()
        except BaseException as err:
            print(f"couldn't disconnect from computer"
                f"\n{type(err).__name__}: {err}")
        try:
            self.cam.disconnect()
        except BaseException as err:
            print(f"couldn't disconnect from Flir camera"
                f"\n{type(err).__name__}: {err}")
        try:
            (self.mog
                .set_frequency(1, 90.0).set_power(1, 29.04)
                .set_frequency(3, f0)
                .set_frequency(4, f0)
                .disconnect()
            )
        except BaseException as err:
            print(f"couldn't reset MOGRF"
                f"\n{type(err).__name__}: {err}")

    def postcmd(self, *args):
        if take_background:
            self.names.append("background")
        arrays = dict(zip(self.names[1:], self.frames[1:]))
        data = FluorescenceScanData(
            outdir=outdir,
            arrays=arrays,
            config=camera_config,
            comments=comments
        )
        data.compute_results(size_fit=False, subtract_bkgd=True, debug=False)
        data.save(arrays=True)
        #data.render_arrays()

    def cmd_visualize(self, *args):
        seq = make_sequence("")
        try:
            tmin, tmax = float(args[0]), float(args[1])
        except IndexError:
            tmin, tmax = 0.0, seq.max_time()
        seq.draw_detailed().set_xlim(tmin, tmax).show().close()
        sys.exit(0)

if __name__ == "__main__":
    NarrowCoolingTweezerAlignment().RUN()
