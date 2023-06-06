from math import *
from PIL import Image


def resize_image(input_image_path, output_image_path, new_width, new_height):
    # Open the image file
    with Image.open(input_image_path) as image:
        # Resize the image
        resized_image = image.resize((new_width, new_height))
        # Save the resized image
        resized_image.save(output_image_path)


for i in range(10, 23):
    width = floor(sqrt(2**i))
    filename = "input0_" + str(width) + ".png"
    resize_image("input0_192x192.png", filename, width, width)