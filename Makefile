CFLAGS=-Wall 
LIBS=-O3 -L/usr/X11R6/lib -lm -lpthread -lX11
EXE=SphSPS

all:
	mkdir -p res
	g++ $(CFLAGS) SphSPS.cpp -o $(EXE) $(LIBS)

test:
	./SphSPS -i ./data/test_img1.jpg -k 1000 -m 0.12 -outm test_img1_labelmap.png -outb test_img1_border.png -c ./data/test_img1_contour.png

test_list:
	./scripts/test_list.sh ./data/list_file.txt ./data/ 1000 0.12
