from picamera2 import Picamera2
import time

picam2 = Picamera2()

config = picam2.create_still_configuration(main = {"size": picam2.sensor_resolution})

config['controls'] = { 'ExposureTime':20000000,
                       'AnalogueGain': 1.0,
                       'AwbEnable': False,
                       'ColourGains': [1.5, 1.2],
                       'AeMeteringMode': 1
                       }

picam2.configure(config)
print(config)
print(config['raw'].get('format'))

picam2.start_preview()
picam2.start()
time.sleep(2)

picam2.capture_file("EL_20s_1AG_1A_4.tiff")
picam2.stop_preview()
picam2.stop()