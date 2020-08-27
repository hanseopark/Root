CXX = g++
CFLAGS = -g -O2

all: Address ValPointer ObjPointer

Address: Address.C
	$(CXX) $(CFLAGS) Address.C -o Address

ValPointer: ValPointer.C
	$(CXX) $(CFLAGS) ValPointer.C -o ValPointer

ObjPointer: ObjPointer.C
	$(CXX) $(CFLAGS) ObjPointer.C -o ObjPointer

clean:
	rm -f Address ValPointer ObjPointer
