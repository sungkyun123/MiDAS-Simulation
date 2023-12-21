CC=g++
CFLAGS= -g -O2 -fpermissive
OBJS = ssdsimul.o algorithm.o model.o ssd_config.o hf.o queue.o
TARGET= ssdsimul
LIBS +=\
       -lpthread\

all: $(TARGET)


$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LIBS)

ssdsimul.o: ssdsimul.cpp ssdsimul.h ssd_config.h model.h
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS)

algorithm.o: algorithm.cpp algorithm.h
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS)

model.o: model.cpp model.h 
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS)

ssd_config.o: ssd_config.cpp ssd_config.h
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS)

queue.o: queue.cpp queue.h
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS)

hf.o: hf.cpp hf.h
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS)

clean: 
	rm -f *.o
	rm -f $(TARGET)
	rm -f midas_*
	rm only_model
