from lib import *
from lib import CONNECTIONS as C
import numpy as np
from itertools import product
import sys
import pathlib
import timeit

timestamp = get_timestamp()
print(timestamp)
outdir = DATADIRS.narrow_cooling_tweezer.joinpath(timestamp)

comments = """
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
reps = 50 # repeat shots for statistics
clk_freq = 10e6 # serial_bits clock frequency; Hz
take_background = True # include a background image
flir_two_shots = False # use the Flir to take two shots
check_556_tweezer = False # use the 556 tweezer instead of the probe beams

# timings (not to scale)
#            |                   tau_all                    |
# start | t0 |           tau_ncool          | tau_twzr_load | tau_probe | end
#            |                              |               |
#        <tau_comp>                         |               |
#            |                          <tau_flir>          |
#            |                              |           <tau_andor>
#            | t_ramp | tau_ramp | tau_hold |
#            |
#        <tau_ref>
t0 = 300e-3 # transfer to green MOT; s
tau_flux_block = -15e-3 # time relative to t0 to stop atom flux; s
tau_blue_overlap = 2e-3 # overlap time of blue beams with green beams relative to t0; s
tau_ncool = 100e-3 # narrow cooling/compression time; s
tau_comp = 46e-3 # start compression ramping relative to beginning of narrow cooling; s
T_COMP = t0 + tau_comp + np.linspace(0.0, 4e-3, 51) # compression ramp times
tau_twzr_load = 10e-3 # time to load into tweezers; s
tau_probe = 40e-3 # probe imaging time; s
tau_all = tau_ncool + tau_twzr_load + tau_probe # main sequence time; s

# camera timings
tau_flir = -5e-3 # Flir camera time relative to end of narrow cooling; s
tau_flir_second = 0e-3
tau_andor = -17e-3 # EMCCD camera time relative to end of tweezer loading; s

# coil settings
SHIMS_FB = np.linspace(1.230, 1.230, 1) # Front/Back shim scans: 1.23; +1.2 for 174
SHIMS_LR = np.linspace(0.120, 0.120, 1) # Left/Right shim scans: 0.13; -0.2 for 174
SHIMS_UD = np.linspace(0.510, 0.510, 1) # Up/Down shim scans:    0.53; +0.4 for 174
B_blue = int(441815) # blue MOT gradient setting
B_green = int(44182 * 1.25) # 174: 1.25 -> 1.1
bits_Bset = int(20) # number of bits in a servo setting
bits_DACset = int(4) # number of bits in a DAC-mode setting
B_COMP = np.linspace(B_green, B_green * 1.8, T_COMP.shape[0]) # compression ramp values

# detunings
DET_MOT = np.linspace(0e-3, 0e-3, 1) # green MOT detuning for tweezer loading; MHz
DET_PROBE = 4.07 + np.linspace(-1.5, 1.5, 5) # probe beam or detuning for imaging; MHz
det_block = 10 # shift the RF on the MOT beams to decouple the fiber; MHz
det_556_tweezer = 4.07 # detuning on the 556 tweezer (probe) path; MHz

## MOGRF TABLE PARAMETERS
# general
ramp_N = 1000 # number of steps in the ramp

# timings
tau_ref = -10e-3 # lead time for mogtable to start up; s
tau_ramp = 30e-3 # ramp duration; s
t_ramp = 20e-3 # start of ramp after t0; s
tau_hold = 50e-3 # post-ramp holding time; s

# frequency parameters
f_ramp = 90.0 # start of ramp; MHz
nu_ramp = 3.4 # extent of ramp; MHz

# power parameters
p_ramp_beg = 29.0 # start of ramp; dBm
p_ramp_end = 0.0 # end of ramp; dBm
p_probe = 17.0 # power in probe beam AOMs; dBm
p_556_tweezer = -16.0 # power in the 556 tweezer (probe) AOM; dBm

def make_sequence(name: str, shim_fb: float, shim_lr: float, shim_ud: float,
        rep: int) \
    -> SuperSequence:
    SEQ = SuperSequence(
        outdir.joinpath("sequences"),
        name,
        {
            "Flir": (Sequence()
                + Sequence.digital_pulse(
                    *C.flir_trig,
                    t0 + tau_ncool + tau_flir,
                    camera_config["exposure_time"] * 1e-6 # convert back us -> s
                )
            ).with_color("C2"),

            "EMCCD": (Sequence()
                + Sequence.digital_pulse(
                    *C.andor_trig,
                    t0 + tau_ncool + tau_twzr_load + tau_andor - 27e-3, # shutter time
                    10e-3 # pulse length is arbitrary; currently set to exposure time
                )
            ).with_color("C2"),

            "Block atom flux": (Sequence()
                + Sequence.digital_lohi(
                    *C.push_aom,
                    t0 + tau_flux_block,
                    t0 + tau_all
                )
                + Sequence.digital_lohi(
                    *C.push_sh,
                    t0 + tau_flux_block - 15e-3, # turn off the shutter ahead of AOM
                    t0 + tau_all
                )
                + Sequence.digital_lohi(
                    *C.mot2_blue_aom, # remove atoms from push beam path to be sure
                    t0 + tau_flux_block,
                    t0 + tau_all
                )
            ).with_color("C3"),

            "CMOT servo ramp": (Sequence.joinall(*[
                Sequence.serial_bits_c(
                    C.mot3_coils_sig,
                    t_comp,
                    int(b_comp), bits_Bset,
                    AD5791_DAC, bits_DACset,
                    C.mot3_coils_clk, C.mot3_coils_sync,
                    clk_freq
                ) for t_comp, b_comp in zip(T_COMP, B_COMP)
            ])).with_color("C1"),

            "MOT servo off": (Sequence()
                + Sequence.serial_bits_c(
                    C.mot3_coils_sig,
                    t0 + tau_ncool, # time for MOT to disperse
                    0, bits_Bset, # turn off to disperse MOT before imaging
                    AD5791_DAC, bits_DACset,
                    C.mot3_coils_clk, C.mot3_coils_sync,
                    clk_freq
                )
            ).with_color("C1"),

            "Blue MOT coil reset": (
                Sequence.serial_bits_c(
                    C.mot3_coils_sig,
                    t0 + tau_all,
                    B_blue, bits_Bset, # reset to normal value for blue MOT
                    AD5791_DAC, bits_DACset,
                    C.mot3_coils_clk, C.mot3_coils_sync,
                    clk_freq
                )
            ).with_color("C1"),

            "Load blue MOT": (Sequence()
                + Sequence.digital_hilo(
                    *C.mot3_blue_aom,
                    0.0, # beams on to load atoms from beginning
                    t0 + tau_blue_overlap 
                )
                + Sequence.digital_hilo(
                    *C.mot3_blue_sh,
                    0.0, # beams on to load atoms from beginning
                    t0 + tau_blue_overlap - 5e-3 # close shutter ahead of AOM
                )
            ).with_color("C0"),

            "FB/LR/UD shims": (Sequence()
                + [
                    Event.analog(**c, s=s) @ (t0 + k * 1e-16) # have to offset times from e/o
                    for k, (c, s) in enumerate([
                        (C.shim_coils_fb, shim_fb),
                        (C.shim_coils_lr, shim_lr),
                        (C.shim_coils_ud, shim_ud)
                    ])
                ]
                + [
                    Event.analog(**c, s=0.0)
                        @ (t0 + tau_ncool + k * 1e-16) # have to offset times from e/o
                    for k, c in enumerate([
                        C.shim_coils_fb, C.shim_coils_lr, C.shim_coils_ud
                    ])
                ]
                + [
                    Event.analog(**c, s=c.default) # have to offset times from e/o
                        @ (t0 + tau_all + k * 1e-16)
                    for k, c in enumerate([
                        C.shim_coils_fb, C.shim_coils_lr, C.shim_coils_ud
                    ])
                ]
            ).with_color("C8"),

            "Green MOT(s)": (Sequence()
                + Sequence.digital_hilo(
                    *C.mot3_green_aom,
                    t0 + tau_ref - 5e-3, # the mogtable is triggered on a falling edge
                    t0 + tau_ref
                )
                + [ # don't use the shutter to turn off the CMOT
                    Event.digital1(
                        **C.mot3_green_sh, s=1)
                    @ (t0 - 5e-3)
                ]
                + Sequence.serial_bits_c(
                    C.mot3_coils_sig,
                    t0,
                    B_green, bits_Bset,
                    AD5791_DAC, bits_DACset,
                    C.mot3_coils_clk, C.mot3_coils_sync,
                    clk_freq
                )
                + Sequence.serial_bits_c(
                    C.mot3_coils_sig,
                    t0 + tau_ncool,
                    0, bits_Bset, # turn off coils for tweezer loading
                    AD5791_DAC, bits_DACset,
                    C.mot3_coils_clk, C.mot3_coils_sync,
                    clk_freq
                )
            ).with_color("C6"),

            "Scope": (
                Sequence.digital_hilo(
                    *C.scope_trig,
                    t0,
                    t0 + tau_ncool
                )
            ).with_color("C7"),

        },
        CONNECTIONS
    )
    SEQ["Sequence"] = ( # dummy sequence; ensures that EW backend stays in-sequence
        Sequence.digital_hilo(
            *C.dummy,
            0.0,
            t0 + tau_all + 1e-3
        )
    ).with_color("k")
    if flir_two_shots:
        SEQ["Flir"] = (
            SEQ["Flir"]
            + Sequence.digital_pulse(
                *C.flir_trig,
                t0 + tau_ncool + tau_flir_second,
                camera_config["exposure_time"] * 1e-6 # convert us -> s
            )
        ).with_color("C2")
    if check_556_tweezer:
        SEQ["556 tweezer"] = (Sequence()
            + Sequence.digital_pulse(
                *C.tweezer_green_aom,
                t0 + tau_ncool,
                camera_config["exposure_time"] * 1e-6 # convert us -> s
                    + 2e-3 # make sure the tweezer is on while Flir is exposed
            )
            # + Sequence.digital_pulse(
            #     *C.tweezer_green_sh,
            #     t0 + tau_ncool - 5e-3, # shutter time and a bit extra for noise
            #     camera_config["exposure_time"] * 1e-6 + 10e-3
            # )
        ).with_color("g")
    else:
        SEQ["Green imaging"] = (Sequence()
            + Sequence.digital_pulse(
                *C.raman_green_aom,
                t0 + tau_ncool,
                tau_twzr_load + tau_probe
            )
            # + Sequence.digital_pulse(
            #     *C.raman_green_sh,
            #     t0 + tau_ncool - 2e-3, # account for 2 ms shutter delay
            #     tau_twzr_load + tau_probe
            # )
        ).with_color("C6")

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

        # "Scope": (Sequence
        #     .digital_hilo(*C.scope_trig, 0.0, 600e-3)
        # ).with_color("C7"),

    }, CONNECTIONS)

# DEFINE MOGRF SEQUENCE
def make_mot_mogtable(mot_det: float):
    # ramp_N is number of steps, not points, so use ramp_N + 1 for number of points
    RAMP_T = np.linspace(t_ramp - tau_ref, t_ramp - tau_ref + tau_ramp, ramp_N + 1)
    RAMP_F = np.linspace(f_ramp, f_ramp + nu_ramp, ramp_N + 1)
    RAMP_P = np.linspace(p_ramp_beg, p_ramp_end, ramp_N + 1)

    mogtable = MOGTable()

    mogtable \
        << MOGEvent(frequency=RAMP_F[0], power=RAMP_P[0]) \
            @ (0.0 - tau_ref) # make sure the frequency/power is right at t0
    for t, f, p in zip(RAMP_T, RAMP_F, RAMP_P):
        mogtable << MOGEvent(frequency=f, power=p) @ t
    # detune and depower the MOT beams for loading tweezers after the CMOT ramp
    #   is finished
    mogtable \
        << MOGEvent(frequency=RAMP_F[-1] + det_block, power=-50.0) \
            @ (RAMP_T[-1] + tau_hold + 1.7e-3) # unknown source of timing error

    return mogtable

# SCRIPT CONTROLLER
class NarrowCoolingTweezerAlignment(Controller):
    def precmd(self, *args):
        self.comp = MAIN.connect()
        self.comp.set_defaults().def_digital(*C.raman_green_sh, 1)
        if check_556_tweezer:
            self.comp.def_digital(*C.tweezer_green_sh, 1)

        self.cam = FLIR.connect()
        (self.cam
            .configure_capture(**camera_config)
            .configure_trigger()
        )

        self.mog = MOGRF.connect()
        if check_556_tweezer:
            self.mog.set_frequency(2, 90.0 + det_556_tweezer).set_power(2, p_556_tweezer)

        # pre-construct all sequences and mogtables so we don't have to do it
        #   in the main loop
        self.names = list()
        self.sequences = list()
        self.ssequences = list()
        mogparams = list(product(DET_MOT, DET_PROBE))
        ewparams = list(product(SHIMS_FB, SHIMS_LR, SHIMS_UD))
        for d_mot, d_probe in mogparams:
            for fb, lr, ud in ewparams:
                for rep in range(reps):
                    name = f"probe-det={d_probe:.5f}_fb={fb:.5f}_lr={lr:.5f}_ud={ud:.5f}_{rep}"
                    self.names.append(name)
        for fb, lr, ud in ewparams:
            for rep in range(reps):
                sseq = make_sequence(name, fb, lr, ud, rep)
                self.sequences.append(sseq.to_sequence())
                self.ssequences.append(sseq)

        # track the mogtables separately because they're slow to load to the controller
        # they should be the slowest-varying parameter in scans
        self.mot_mogtables = [make_mot_mogtable(d_mot) for d_mot in DET_MOT]

    def run_sequence(self, *args):
        for table in self.mot_mogtables:
            (self.mog
                .set_frequency(1, 90.0).set_power(1, -50.0) # make sure the AOM is off
                .set_mode(1, "TSB")
                .table_load(1, table)
                .table_arm(1)
                .set_table_rearm(1, True)
            )
            for d_probe in DET_PROBE:
                (self.mog
                    .set_frequency(3, 90.0 + d_probe).set_power(3, p_probe)
                    .set_frequency(4, 90.0 + d_probe).set_power(4, p_probe)
                )
                for seq in self.sequences:
                    (self.comp
                        .enqueue(seq)
                        .run()
                        .clear()
                    )
            (self.mog
                .table_stop(1)
                .table_clear(1)
            )
        if take_background:
            (self.comp
                .enqueue(seq_bkgd.to_sequence())
                .run()
                .clear()
            )

        self.comp.clear().disconnect()

        (self.mog
            .table_stop(1)
            .table_clear(1)
            .set_frequency(1, 90.0).set_power(1, 29.04)
            .set_frequency(3, 94.07)
            .set_frequency(4, 94.07)
            .set_mode(1, "NSB")
		    .set_output(1, True)
            .disconnect()
        )

    def run_camera(self, *args):
        self.frames = self.cam.acquire_frames(
            num_frames = (
                len(self.names)
                + int(flir_two_shots) * len(self.names)
                + int(take_background)
            ),
            timeout=5, # s
            roi=[968, 737, 40, 40]
        )
        # might be important to disconnect in the action rather than postcmd
        #   due to possible bad sharing of resources between threads
        self.cam.disconnect()

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
                .table_stop(1)
                .table_clear(1)
                .set_frequency(1, 90.0).set_power(1, 29.04)
                .set_frequency(3, 93.98)
                .set_frequency(4, 93.98)
                .set_mode(1, "NSB")
                .set_output(1, True)
                .disconnect()
            )
        except BaseException as err:
            print(f"couldn't reset MOGRF"
                f"\n{type(err).__name__}: {err}")

    def postcmd(self, *args):
        names = list()
        for name in self.names:
            if flir_two_shots:
                avgimgdir = outdir.joinpath("images").joinpath("averages")
                avgimgdir_pre = avgimgdir.joinpath("pre")
                avgimgdir_pre.mkdir(parents=True, exist_ok=True)
                avgimgdir_post = avgimgdir.joinpath("post")
                avgimgdir_post.mkdir(parents=True, exist_ok=True)
                name_split = name.split("_")
                names.append(str(
                    avgimgdir_pre.joinpath(
                        "_".join(name_split[:-1]) + "_pre_" + name_split[-1])
                ))
                names.append(str(
                    avgimgdir_post.joinpath(
                        "_".join(name_split[:-1]) + "_post_" + name_split[-1])
                ))
            else:
                names.append(name)
        if take_background:
            names.append("background")
        arrays = {
            name: frame for name, frame in zip(names, self.frames)
        }
        data = FluorescenceScanData(
            outdir=outdir,
            arrays=arrays,
            config=camera_config,
            comments=comments
        )
        data.compute_results(
            size_fit=False,
            subtract_bkgd=True,
            debug=False,
            mot_number_params={
                # from measurement on 02.21.22
                "intensity_parameter": 2 * 3.55,
                # from measurement on 02.15.22
                "detuning": abs(93.4e6 - 94.08e6) * 2 * np.pi,
            }
        )
        data.save(arrays=True)
        #data.render_arrays()

    def cmd_visualize(self, *args):
        seq = make_sequence("", SHIMS_FB.max(), SHIMS_LR.max(), SHIMS_UD.max(), 0)
        try:
            tmin, tmax = float(args[0]), float(args[1])
        except IndexError:
            tmin, tmax = 0.0, seq.max_time()
        P = seq.draw_detailed(mogtables=[
            ( make_mot_mogtable(0.0).with_color("C6"),
                dict(name="Green MOT AOM", offset=t0 + tau_ref) ),
        ])
        for t in [
            t0,
            t0 + tau_ncool,
            t0 + tau_ncool + tau_twzr_load,
            t0 + tau_ncool + tau_twzr_load + tau_probe,
        ]:
            P.axvline(t, color="r", linestyle="-", linewidth=0.4)
        for t in [
            t0 + t_ramp,
            t0 + t_ramp + tau_ramp,
            t0 + t_ramp + tau_ramp + tau_hold,
        ]:
            P.axvline(t, color="g", linestyle="--", linewidth=0.4)
        for t in [
            t0 + tau_comp,
            t0 + tau_ncool + tau_flir,
            t0 + tau_ncool + tau_twzr_load + tau_andor,
        ]:
            P.axvline(t, color="b", linestyle=":", linewidth=0.4)
        (P
            .set_xlim(tmin, tmax)
            .show()
            .close()
        )
        sys.exit(0)

if __name__ == "__main__":
    NarrowCoolingTweezerAlignment().RUN()
