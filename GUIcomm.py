import pickle
from  lib import *
from PIL import Image

import cmot_tof, freespace_green_resonance, narrow_cooling_tweezer_align

'''
"cmot_tof.py" =="CMTF"
"freespace_green_resonance.py" =="FSGR"
"narrow_cooling_tweezer_align.py" == "NCTA"


'''
outdir=""

def create_ss_obj(tag, params = 0.0):
    if (tag == "CMTF"):
        SEQ = cmot_tof.make_sequence("seq", params)
        return SEQ
        
    if (tag == "FSGR"):
        SEQ = freespace_green_resonance.make_sequence("seq")
        return SEQ

    if (tag == "NCTA"):
        SEQ = narrow_cooling_tweezer_align.make_sequence("seq", params)
        return SEQ


def sequences_in_super(Supersequence):
    sequences_dict = Supersequence.to_primitives()
    return sequences_dict

def save_altered_super(seqs_dict):
    new_super_seq = SuperSequence.from_primitives(outdir, "new_SSS", seqs_dict)
    new_super_seq.save()

def events_in_seq(name, seq):

    supsequence = SuperSequence(outdir, name, seq)
    sequenc = supsequence.by_times()
    print(sequenc)
    print(type(sequenc))
    
    '''
    sequence = ssequence.to_sequence()
    events_dict = sequence.to_primitives()
    for i in events_dict.keys():
        print(i)
    return events_dict
    '''

 




def save_object(obj, filename):
    with open(filename, 'wb') as outp:
        pickle.dump(obj, outp, pickle.HIGHEST_PROTOCOL)


def create_dig_event(board, channel, voltage, state, name, time):
    eventaddr = 'Events/' + name + '_dig' + '.txt'
    newevent = open(eventaddr, 'w')

    #Fix string -> ints for mask and state
    event = Event(EventType.Digital, board, connector=0, mask=ch(), state=ch())
    save_object(event, eventaddr)
    



def create_ana_event(board, channel, voltage, state, name, time):
    eventaddr = 'Events/' + name + '_ana' + '.txt'
    newevent = open(eventaddr, 'w')


    event = Event(EventType.Analog, board, connector=0, mask=ch(), state=ch())
    save_object(event, eventaddr)


def create_flir_seq():
    pass

def create_EMCCD_seq():
    pass

def create_atom_flux_seq():
    pass

def create_CMOT_servo_seq():
    pass

def create_bMOT_coil_seq():
    pass

def create_bMOT_load_seq():
    pass

def create_gMOT_seq():
    pass

def create_scope_seq():
    pass


def create_supersequence():
    pass
