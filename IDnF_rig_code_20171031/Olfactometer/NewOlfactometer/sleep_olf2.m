function [] = sleep_olf2()

close_serial_port(13)
olf_arduino_serial_comm(0, [.5, .5], 1, 1)
pause(0.5)
trigger_scan(1)
pause(0.5)
trigger_scan(0)
olf_arduino_serial_comm(2, [.5, .5], 1, 1)
pause(0.5)
close_serial_port(13)
trigger_scan(1)
pause(0.5)
trigger_scan(0)