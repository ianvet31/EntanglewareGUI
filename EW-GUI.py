from tkinter import HORIZONTAL
from lib.ew.lib.entangleware_control_link import set_analog_state
from lib.rtcontrol.driver import set_digital_state
from lib.rtcontrol.parser import CONNECTOR_D_RANGE
from lib.system import MAIN
from lib.ew import EventType, DigitalConnection, AnalogConnection
from lib.rtcontrol import *

from email.policy import default
from turtle import update
from venv import create
import dearpygui.dearpygui as dpg
import matplotlib.pyplot as plt

import numpy as np
from GUIcomm import *


dpg.create_context()
dpg.setup_dearpygui()

dpg.enable_docking(dock_space=True)
dpg.create_viewport(title='Entangleware Control GUI', width=1600, height=900)
dpg.show_viewport()

counter = 0


class Sequences():
    sequences = []



"""

texture_data = []
for i in range(0, 100 * 100):
    texture_data.append(255 / 255)
    texture_data.append(0)
    texture_data.append(255 / 255)
    texture_data.append(255 / 255)

with dpg.texture_registry(show=True):
    dpg.add_dynamic_texture(100, 100, texture_data, tag="texture_tag")


def _update_dynamic_textures(sender, app_data, user_data):
    new_color = dpg.get_value(sender)
    new_color[0] = new_color[0] / 255
    new_color[1] = new_color[1] / 255
    new_color[2] = new_color[2] / 255
    new_color[3] = new_color[3] / 255

    new_texture_data = []
    for i in range(0, 100 * 100):
        new_texture_data.append(new_color[0])
        new_texture_data.append(new_color[1])
        new_texture_data.append(new_color[2])
        new_texture_data.append(new_color[3])

    dpg.set_value("texture_tag", new_texture_data)


# for inside window,
    dpg.add_image("texture_tag")
    dpg.add_color_picker((255, 0, 255, 255), label="Texture",
                         no_side_preview=True, alpha_bar=True, width=200,
                         callback=_update_dynamic_textures)


"""

#Eventbuilder vars

signal_type = 0    #digital = 0, analog = 1
voltage = 0.0
board = 0
channel = 0
time = 0.0
state = ""
name = ""

#Eventbuilder callbacks
def CBsignal_type(sender, data):
    if data == "Analog":
        signal_type = 1
    else:
        signal_type = 0

    print(signal_type)
        

def CBvoltage(sender, data):
    voltage = data

def CBboard(sender, data):
    board = data

def CBtime(sender, data):
    time = data

def CBstate(sender, data):
    state = data

def CBname(sender, data):
    name = data

def CBsave_event(sender, data):
    if (signal_type == 1):
        event = create_ana_event(board, channel, voltage, state, name, time)
    else:
        event = create_dig_event(board, channel, voltage, state, name, time)

#Cam callbacks

def CBcam_settings():
    pass
def CBExposure_Time():
    pass
def CBGain():
    pass
def CBGamma():
    pass
def CBBlack_Level():
    pass
def CBBin_Size():
    pass

#Sequencevars

currentevents = {}

#Sequencebuilder callbacks
def load_sequences(sequences):
    pass


def CBsave_sequence(sender, data):
    pass



#Supersequencebuilder callbacks
def CBpreset_seq(sender, data):
    pass



'''

width, height, channels, data = dpg.load_image("Plot.png")
with dpg.texture_registry(show=False):
    dpg.add_dynamic_texture(width, height, data, tag="plot_tg")

with dpg.window(label="Sequence Plot", tag="Plotter", pos=[0, 400], width=1200, height=500):
        dpg.add_image("plot_tg")

'''


def CB_render_supersequence(sender, value, user_data):
    """
    render_sequence(user_data)
    width, height, channels, data = dpg.load_image("Plot.png")
    dpg.set_value("plot_tg", data)
    """
    pass


def event_saved(sender, user_data):

    _e_type = dpg.get_value("event_signal")
    if (_e_type == "Digital"):
        _e_type = True
    else:
        _e_type = False
    
    _e_time = dpg.get_value("event_time")
    _e_conn = dpg.get_value("event_conn")
    _e_mask = dpg.get_value("event_mask")
    _e_state = dpg.get_value("event_state")
    _e_board = dpg.get_value("event_board")

    
    #Only handling digital type events?

    new_event = Event(EventType.Digital, time=_e_time, board=_e_board, connector=_e_conn, mask=_e_mask, state=_e_state)


    # retrives dictionary indices of current seq and event to save the new event to the correct spot


    event_idx = dpg.get_item_user_data("event_save")
    

    seqs_dict_old = dpg.get_item_user_data("seqs_listbox")
    seq_old = seqs_dict_old[dpg.get_item_user_data("seq_save")]
    events_lst = seq_old["events"]
    
    new_events = events_lst
    new_events[event_idx] = new_event.to_primitives()
    
    seq_new = seq_old
    seq_old["events"] = seq_new
    seqs_dict_old[dpg.get_item_user_data("seq_save")] = seq_old


    pass


def E_clicked(sender, value, user_data):

    # retrieves index to display order and saves for future use

    string = dpg.get_value(sender)
    idx = string[6]
    idx = int(idx)
    dpg.configure_item("event_num", default_value = "Event " + str(idx))
    dpg.set_item_user_data("event_save", user_data=idx-1)
    


    # extracts event's values for event builder (implement for analog!)

    event = user_data[idx - 1]
    e_type = event["is_digital"]
    if (e_type == True):
        e_type = "Digital"
    else:
        e_type = "Analog"


    e_time = event["time"]
    e_conn = event["connector"]
    e_mask = event["mask"]
    e_state = event["state"]
    if ("board" in event):
            e_board = event["board"]
            dpg.configure_item("event_board", default_value = e_board)

    dpg.configure_item("event_time", default_value = e_time)
    dpg.configure_item("event_signal", default_value = e_type)
    dpg.configure_item("event_mask", default_value = e_mask)
    dpg.configure_item("event_conn", default_value = e_conn)
    dpg.configure_item("event_state", default_value = e_state)

    
def S_clicked(sender, value, user_data):


    dpg.set_item_user_data("seq_save", dpg.get_value(sender))
    


    curr_seq = user_data[dpg.get_value(sender)]
    events_lst = curr_seq["events"]

    items=[]
    counter=1
    for i in events_lst:
        e_time = i["time"]
        items.append("Event " + str(counter) + ", at time " + str(e_time))
        counter += 1

    dpg.configure_item("events_listbox", items=items)
    dpg.set_item_user_data("events_listbox", events_lst)
    
    


def SS_clicked(sender, value):

    SS = create_ss_obj(sender)

    seqs_dct = sequences_in_super(SS)
    
    seqs_keys = list(seqs_dct.keys())

    dpg.configure_item("seqs_listbox", items=seqs_keys)
    dpg.set_item_user_data("seqs_listbox", seqs_dct)

    dpg.set_item_user_data("ss_run", SS)



# Decide on useful parameters for CMTF supersequence, some constants won't be useful to adjust during testing...

def cmtf_params_clicked():
    with dpg.window(label = "CMTF Config", pos=[100, 50], width=300):
        dpg.add_text("Camera Config")
        dpg.add_input_float(label="Exposure Time", default_value = 300.0, callback=CBExposure_Time)
        dpg.add_input_float(label="Gain", default_value=47.99, callback=CBGain)
        dpg.add_input_float(label="Gamma", default_value = 1.25, callback=CBGamma)
        dpg.add_input_float(label="Black Level", default_value = 5.0, callback=CBBlack_Level)
        dpg.add_input_int(label="Bin Size", default_value = 1,  callback=CBBin_Size)

        dpg.add_spacing(count=5)
        dpg.add_text("Main Sequence Parameters")
        dpg.add_input_int(label="Repeated shots", default_value=5)

        dpg.add_input_float(label="Transfer to Green MOT (t0)", default_value = 300e-3)
        dpg.add_input_float(label="Tau Flux Block", default_value=-15e-3)
        dpg.add_input_float(label="Tau Blue Overlap", default_value=2e-3)
        dpg.add_input_float(label="Cooling/Compression time", default_value= 100e-3)
        dpg.add_input_float(label="Tau compression start", default_value= 46e-3)
        dpg.add_input_float(label="Probe imaging time", default_value=0.5e-3)

        dpg.add_spacing(count=5)
        dpg.add_text("Coil Settings")

        dpg.add_input_float(label="Shim FB", default_value=1.275)
        dpg.add_input_float(label="Shim LR", default_value=0.205)
        dpg.add_input_float(label="Shim UD", default_value=0.528)
        dpg.add_input_int(label="B Blue", default_value=441815)
        dpg.add_input_int(label="B Green", default_value = int(44182*1.25))

        dpg.add_button(label="Save")



# Decide on useful parameters for FSGR supersequence, some constants won't be useful to adjust during testing...


def fsgr_params_clicked():
    with dpg.window(label = "CMTF Config", pos=[100, 50]):
        dpg.add_text("Camera Config")
        dpg.add_input_float(label="Exposure Time", default_value = 300.0, callback=CBExposure_Time)
        dpg.add_input_float(label="Gain", default_value=47.99, callback=CBGain)
        dpg.add_input_float(label="Gamma", default_value = 1.25, callback=CBGamma)
        dpg.add_input_float(label="Black Level", default_value = 5.0, callback=CBBlack_Level)
        dpg.add_input_int(label="Bin Size", default_value = 1,  callback=CBBin_Size)
        
        dpg.add_spacing(count=5)
        dpg.add_text("Main Sequence Parameters")
        dpg.add_input_int(label="Repeated shots", default_value=15)

        dpg.add_input_float(label="Release Blue MOT (t0)", default_value = 450e-3)
        dpg.add_input_float(label="Tau Flux Block", default_value=-15e-3)
        dpg.add_input_float(label="Tau TOF", default_value=0e-3)
        dpg.add_input_float(label="Probe imaging after release", default_value=0.5e-3)

        dpg.add_spacing(count=5)
        dpg.add_text("Coil Settings")

        dpg.add_input_int(label="B blue", default_value=441815)
        dpg.add_input_int(label="# of bits in servo", default_value=20)
        dpg.add_input_int(label="# of bits in DAC", default_value=4)
        dpg.add_button(label="Save")
 

def _config(sender, keyword, user_data):

    widget_type = dpg.get_item_type(sender)
    items = user_data

    if widget_type == "mvAppItemType::mvRadioButton":
        value = True

    else:
        keyword = dpg.get_item_label(sender)
        value = dpg.get_value(sender)

    if isinstance(user_data, list):
        for item in items:
            dpg.configure_item(item, **{keyword: value})
    else:
        dpg.configure_item(items, **{keyword: value})



def RT_callback():
    rttype = dpg.get_value("RT_type")
    rtconn = dpg.get_value("RT_conn")
    rtmask = dpg.get_value("RT_mask")
    rtstate = dpg.get_value("RT_state")

    if (rttype == "Digital"):
        set_digital_state([rtconn, rtmask, int(rtstate)])
    if (rttype == "Analog"):
        set_analog_state([rtconn, rtstate])
    else:
        pass


with dpg.window(label="Supersequences", pos=[0,0], width=1600, height = 400, no_resize=True, no_move=True):
    with dpg.group():
        dpg.add_spacing(count = 5)
        dpg.add_text("Preset Super Sequences")
        with dpg.group(horizontal=True, horizontal_spacing=170):

            dpg.add_button(label="CMOT_TOF", tag="CMTF", callback=SS_clicked) 
            dpg.add_button(label="CMTF parameters", tag="CMTFp",callback=cmtf_params_clicked)

        with dpg.group(horizontal=True, horizontal_spacing=53):

            dpg.add_button(label="Freespace Green Resonance", tag="FSGR", callback=SS_clicked)
            dpg.add_button(label="FSGR parameters", tag="FSGRp", callback=fsgr_params_clicked)

       
        
        dpg.add_spacing(count = 5)

        dpg.add_button(label="Save")
        dpg.add_button(label="Load")
        dpg.add_button(label="Run", tag="ss_run")

    with dpg.child_window(pos = [400, 0], width = 250, height = 400):
        dpg.add_spacing(count=5)
        dpg.add_text("Sequences")
        dpg.add_spacing(count=5)
        
        dpg.add_listbox(tag="seqs_listbox", num_items = 15, width = 250, callback=S_clicked)
        

    with dpg.child_window(pos=[650,0], width=300, height=400):

        dpg.add_spacing(count=5)
        dpg.add_text("Events")
        dpg.add_spacing(count=5)

        dpg.add_listbox(tag="events_listbox", width = 300, num_items = 15, callback=E_clicked)
        dpg.add_spacing(count=5)
        dpg.add_button(label="Save", tag="seq_save")

    with dpg.child_window(pos=[950, 0], width=550, height = 400):
        
        dpg.add_spacing(count=5)
        dpg.add_text("Event Builder")

        dpg.add_spacing(count=5)
        dpg.add_text("", tag="event_num")
        
        dpg.add_radio_button(items=["Digital", "Analog"], label="Signal Type", tag="event_signal", default_value=None)
        dpg.add_spacing(count = 5)

        dpg.add_input_float(label="Time", tag="event_time")

       
        dpg.add_input_int(label = "Board", tag="event_board")
    
        dpg.add_input_int(label = "Connector", tag="event_conn")
        dpg.add_spacing(count = 1)
        dpg.add_input_int(label="Mask", tag="event_mask")
        dpg.add_spacing(count=1)
        dpg.add_input_int(label="State", tag="event_state")

        dpg.add_spacing(count = 5)
        dpg.add_button(label = "Save", tag="event_save", callback=event_saved)
        

with dpg.window(label="Real time control", pos=[800, 0], width=400):
    dpg.add_radio_button(items=["Digital", "Analog"], label="Type", tag="RT_type")
    dpg.add_spacing(count=5)
    dpg.add_input_int(label="Connector", max_value=4, tag="RT_conn")
    dpg.add_input_int(label="Mask (Only for digital)", tag="RT_mask")
    dpg.add_input_float(label="State", tag="RT_state")
    dpg.add_button(label="Run", callback=RT_callback)
    



#Plot Rendering

with dpg.window(label="Sequence Plot", pos=[0, 400], width=1600, tag="sequence_plot"):
    with dpg.plot(label="Sequence Plot",width=1600, height=400):

        SS = create_ss_obj(tag = "CMTF")

        # SS = dpg.get_item_user_data("ss_run")
        
        _connections = list(SS.defaults.connections.items())

        T0 = min([seq.min_time() for seq in SS.values()])
        T1 = max([seq.max_time() for seq in SS.values()])
        color_C = 0
        hk = 3
        linewidth=0.7

        def ycoord(x, k):
            y = x
            return hk * k + 1 + y


        conns_list=[]
        xs=[]
        xs2=[]

        xs_main = []
       
        for k, (conn_label, conn) in enumerate(reversed(_connections)):
                conns_list.append((conn_label, ycoord(0.5, k)))
                xs.append(ycoord(0, k))
                xs2.append(ycoord(1, k))
                xs_main.append(ycoord(0.5, k))

        conns_list = tuple(conns_list)

        dpg.add_plot_legend()
        dpg.add_plot_axis(dpg.mvXAxis, label="time", tag="time_axis")
        dpg.add_plot_axis(dpg.mvYAxis, label="channel", tag="channel_axis")
        dpg.set_axis_ticks(dpg.last_item(), conns_list)


        '''
        dpg.add_hline_series(xs, parent="channel_axis")
        dpg.add_hline_series(xs2, parent="channel_axis")
        '''
        dpg.add_hline_series(xs_main, parent = "channel_axis")
        
        

        for seq_label, seq in SS.by_times_named():
            if seq.color is None:
                color = f"C{color_C}"
                color_C = (color_C + 1) % 10
            else:
                color = seq.color
            t0 = seq.min_time()
            dpg.add_plot_annotation(label=seq_label, default_value = (t0, 0), color=[0, 0, 0, 0])
            dpg.add_vline_series([t0], parent="channel_axis")
            for k, (conn_label, conn) in enumerate(reversed(_connections)):
                ts = seq._get_states(**conn)
                if len(ts) > 0:
                    dpg.add_plot_annotation(label = seq_label, default_value=(t0, hk * k + 1), color=[0, 0, 0, 0])


    
    
    


# below loop is interchangeable with dpg.start_dearpygui()

#while dpg.is_dearpygui_running():
#    dpg.render_dearpygui_frame()

dpg.start_dearpygui()
dpg.destroy_context()
