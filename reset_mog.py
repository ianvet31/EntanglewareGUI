import lib
mog = lib.MOGRF.connect()
(mog
    .table_stop(1)
    .table_clear(1)
    .set_frequency(1, 90.0).set_power(1, 29.04)
    .set_mode(1, "NSB")
    .set_output(1, True)
    .disconnect()
)
