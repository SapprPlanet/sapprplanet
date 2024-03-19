import os
rootdir = '../data/maven.mag.calibrated'

for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        if file.endswith('.sts'):
            filename = os.path.join(subdir, file)
            print(filename)
