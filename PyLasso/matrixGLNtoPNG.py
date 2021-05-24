from PIL import Image, ImageDraw, ImageFont

import sys

filename=sys.argv[1]
#filename = 'proteins/LASSOLINKINGS_1zd0_48_131_A.xyz'

f = open(filename)
line = f.readline() #1.linia
words_1line = line.split()

CHAIN = int(words_1line[1])
#KLATKA = 20
#WIDTH = CHAIN*KLATKA
#HEIGHT = CHAIN*KLATKA

KLATKA = int(1100/CHAIN)
FLOATKLATKA = float(1100)/CHAIN
WIDTH = 1100
HEIGHT = 1100

FONT_SIZE=45
FONT_SIZE_SMALL=20
FONT_SIZE_MIDDLE=30

image = Image.new('RGB', (WIDTH, HEIGHT), 'white')
draw = ImageDraw.Draw(image)

data = f.readlines()

for line in data[1:]:
  words = line.split()
#  print words
#  color = (255,int(255-255*float(words[3])),int(255-255*float(words[3]))) if words[2][0] == "m" else (int(255-255*float(words[3])),int(255-255*float(words[3])),255)
  color = (int(255*float(words[5])),int(255*float(words[6])),int(255*float(words[7])))
  for x in range(0,KLATKA+1):
    for y in range(0,KLATKA+1):
      draw.point((int(int(words[1])*FLOATKLATKA)+x, int(int(words[0])*FLOATKLATKA)+y), fill=color)

font = ImageFont.truetype("arial.ttf", FONT_SIZE)
font_middle = ImageFont.truetype("arial.ttf", FONT_SIZE_MIDDLE)
font_small = ImageFont.truetype("arial.ttf", FONT_SIZE_SMALL)
#font = ImageFont.truetype("/home/shared/baza_lassa/surfacesMyTraj_list2016/arial.ttf", FONT_SIZE)
#font_middle = ImageFont.truetype("/home/shared/baza_lassa/surfacesMyTraj_list2016/arial.ttf", FONT_SIZE_MIDDLE)
#font_small = ImageFont.truetype("/home/shared/baza_lassa/surfacesMyTraj_list2016/arial.ttf", FONT_SIZE_SMALL)

# min i max
max = "%.2f" % float(words_1line[8])
min = "%.2f" % float(words_1line[5])
xmax = int(int(words_1line[7])*FLOATKLATKA) if int(int(words_1line[7])*FLOATKLATKA)<HEIGHT-FONT_SIZE else HEIGHT-FONT_SIZE-5
ymax = int(int(words_1line[6])*FLOATKLATKA) if int(int(words_1line[6])*FLOATKLATKA)<WIDTH-FONT_SIZE else WIDTH-FONT_SIZE-5
xmin = int(int(words_1line[4])*FLOATKLATKA) if int(int(words_1line[4])*FLOATKLATKA)<HEIGHT-FONT_SIZE else HEIGHT-FONT_SIZE-5
ymin = int(int(words_1line[3])*FLOATKLATKA) if int(int(words_1line[3])*FLOATKLATKA)<WIDTH-FONT_SIZE else WIDTH-FONT_SIZE-5
#draw.text((xmax, ymax), max+" ("+words_1line[6]+","+words_1line[7]+")", font=font, fill=(0,0,0,255))
#draw.text((xmin, ymin), min+" ("+words_1line[3]+","+words_1line[4]+")", font=font, fill=(0,0,0,255))
draw.ellipse((xmax-12, ymax-12, xmax+12, ymax+12), fill = 'blue', outline ='black')
draw.ellipse((xmin-12, ymin-12, xmin+12, ymin+12), fill = 'red', outline ='black')
draw.text((.65*WIDTH, HEIGHT/10.), "max GLN = "+max+" ("+words_1line[7]+","+words_1line[6]+")", font=font_middle, fill=(0,0,0,255))
draw.text((.65*WIDTH, 1.5*HEIGHT/10.), "min GLN = "+min+" ("+words_1line[4]+","+words_1line[3]+")", font=font_middle, fill=(0,0,0,255))

draw.text((10,10), filename+", length "+str(CHAIN), font=font, fill=(0,255,0,255))

#osie
step = int(CHAIN/5)/10*10 if CHAIN>150 else int(CHAIN/5)/5*5 if CHAIN>50 else int(CHAIN/5)
for x in range(1,6):
  draw.text((x*step*FLOATKLATKA, HEIGHT-FONT_SIZE_SMALL-2), str(x*step), font=font_small, fill=(50,50,50,255))
  draw.text((0,x*step*FLOATKLATKA), "-"+str(x*step), font=font_small, fill=(50,50,50,255))

image.save(filename+'.png')

image.save(filename+'.png')

f.close()



#for x, y in [(x, y) for x in range(WIDTH) for y in range(HEIGHT)]:
#  color = (136, 54, 240) if (x / 20 + y / 20) % 2 == 0 else (10, 204, 127)
#  draw.point((x, y), fill=color)

#image.save('image.png')




