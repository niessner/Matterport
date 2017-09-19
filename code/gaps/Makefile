# 
# Makefile for GAPS
#



#
# Make targets
#

opt:
	$(MAKE) target "TARGET=$@"

debug:
	$(MAKE) target "TARGET=$@"

mesa:
	$(MAKE) target "TARGET=$@"

clean:
	$(MAKE) target "TARGET=$@"

release: pkgs apps docs makefiles vc
	mkdir -p release
	cp README.txt release
	cp LICENSE.txt release
	cp Makefile release
	cp -r makefiles release
	$(MAKE) target "TARGET=$@"

target: 
	cd pkgs; $(MAKE) $(TARGET)
	cd apps; $(MAKE) $(TARGET)
	# cd vc; $(MAKE) $(TARGET)
	# cd docs; $(MAKE) $(TARGET)







