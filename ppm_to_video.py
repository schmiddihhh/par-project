import cv2
import os

image_folder = '/home/schmiddihhh/Uni/PAR/project/test_images/'
video_name = '/home/schmiddihhh/Uni/PAR/project/test_images/video.avi'

images = [img for img in os.listdir(image_folder) if img.endswith(".ppm")]
frame = cv2.imread(os.path.join(image_folder, images[0]))
height, width, layers = frame.shape
cv2.imshow("window", frame)

video = cv2.VideoWriter(video_name, 0, 1, (width,height))

for image in images:
    video.write(cv2.imread(os.path.join(image_folder, image)))

cv2.destroyAllWindows()
video.release()