from math import *
from PIL import Image
import numpy as np

def resize_image(input_image_path, output_image_path, new_width, new_height):
    # Open the image file
    with Image.open(input_image_path) as image:
        # Resize the image
        resized_image = image.resize((new_width, new_height))
        # Save the resized image
        resized_image.save(output_image_path)


for i in range(10):
    j = 7*(i/9)+10
    width = ceil(sqrt(2**j))
    filename = "input0_" + str(width) + ".png"
    resize_image("input0_192x192.png", filename, width, width)