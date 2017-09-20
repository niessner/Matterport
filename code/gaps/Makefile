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

target: 
	cd pkgs; $(MAKE) $(TARGET)
	cd apps; $(MAKE) $(TARGET)







