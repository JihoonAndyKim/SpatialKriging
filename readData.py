import os
for filename in os.listdir(os.getcwd()+"/obs"):
    with open("obs/" + filename) as f:
        lines = f.readlines()
    print lines
