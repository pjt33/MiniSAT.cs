EXE = minisat
DLL = MiniSatCS
DLL_SRC = solver.cs structs.cs
CSC = gmcs
DOTNET_OPT = mono --optimize=all
DOTNET_DEBUG = mono --debug

all: $(EXE).exe

$(EXE).exe: $(DLL).dll main.cs
	$(CSC) /debug+ /define:DEBUG main.cs /r:$(DLL).dll /out:$@

$(DLL).dll: $(DLL_SRC)
	$(CSC) /t:library /debug+ /define:DEBUG $(DLL_SRC) /out:$@

$(EXE)-opt.exe: main.cs
	$(CSC) /optimize+ main.cs $(DLL_SRC) /out:$@

bench: $(EXE)-opt.exe
	for n in 1 2 3 4 5 ; do \
	$(DOTNET_OPT) $(EXE)-opt.exe tests/flat200-59.cnf | grep CPU ; \
	done

t: $(EXE).exe
	@$(DOTNET_DEBUG) $(EXE).exe -s tests/sat/*
	@echo " SAT OK"
	@mono $(EXE).exe -u tests/unsat/*
	@echo " UNSAT OK"

clean:
	rm -f *.exe *.mdb *.pdb *.dll
