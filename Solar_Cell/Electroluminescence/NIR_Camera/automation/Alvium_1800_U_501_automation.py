from Aluvium_1800_U_lib import *

## input to the program (folder name, file name, save type)
## camera setting (Exposure, gain)
exposure = 0.2e6 #us
gain = 0 #dB
frame_rate = 5 #fps
folder = "C:\\Drive\\semester 7\\BTP-1\\work4_EL\\Noir_Camera\\automation\\img_temp"
filename = f"\\img_exp_{exposure}_gain_{gain}_fps_{frame_rate}"
savetype = ".png"
img_path = folder + filename + savetype

def main():
    cam_id = parse_args()

    with VmbSystem.get_instance():
        with get_camera(cam_id) as cam:
            # setup general camera settings and the pixel format in which frames are recorded
            setup_camera(cam, exposure_us = exposure ,gain_db = gain, frame_rate = frame_rate)
            setup_pixel_format(cam)
            handler = Handler()

            try:
                # Start Streaming with a custom a buffer of 10 Frames (defaults to 5)
                cam.start_streaming(handler=handler, buffer_count=10)
                print(img_path)
                display = handler.get_image()
                cv2.imwrite(img_path, display)
                print(f"Image saved as {filename + savetype}")
                cam.stop_streaming()
            except:
                print("failed to capture image")
                
if __name__ == '__main__':
    main()
