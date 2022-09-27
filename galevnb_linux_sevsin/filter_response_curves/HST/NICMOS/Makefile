.PHONY: delete


nothing:


all: download fits2data response final rename calib delete 


delete:
	rm *.fits *.data response* synphot1.tar.gz filter.raw

download:
	wget http://stsdas.stsci.edu/synphot/synphot1.tar.gz
	tar xf synphot1.tar.gz
	mv comp/nicmos/*.fits .
	rm -rf comp mtab


conv = /work/galev/filter/ConvolveQE/ConvolveQE
define convolve
	$(conv) $1 $2 $1.$3
endef




define convolve_nic1
	$(call convolve,nic1_f090m_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f095n_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f097n_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f108n_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f110m_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f110w_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f113n_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f140w_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f145m_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f160w_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f164n_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f165m_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f166n_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f170m_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f187n_012_syn.fits.data,response_nic1.res,$1)
	$(call convolve,nic1_f190n_012_syn.fits.data,response_nic1.res,$1)
endef

define convolve_nic2
	$(call convolve,nic2_f110w_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f160w_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f165m_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f171m_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f180m_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f187n_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f187w_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f190n_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f204m_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f205w_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f207m_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f212n_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f215n_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f216n_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f222m_012_syn.fits.data,response_nic2.res,$1)
	$(call convolve,nic2_f237m_012_syn.fits.data,response_nic2.res,$1)
endef

define convolve_nic3
	$(call convolve,nic3_f108n_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f110w_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f113n_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f150w_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f160w_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f164n_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f166n_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f175w_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f187n_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f190n_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f196n_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f200n_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f205m_005_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f212n_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f215n_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f222m_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_f240m_012_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_g096_005_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_g141_005_syn.fits.data,response_nic3.res,$1)
	$(call convolve,nic3_g206_005_syn.fits.data,response_nic3.res,$1)
endef

fits2data:
	/work/python/read_table.py data *.fits


response:
	/work/galev/spec_tools/make_lambdalist 1000 30000 5 | awk '{ print $$1,1; }' >response_stage1
	$(conv) response_stage1 nic_primary_001_syn.fits.data  response_stage2
	$(conv) response_stage2 nic_secondary_001_syn.fits.data response_stage3
	$(conv) response_stage3 nic_reimag_004_syn.fits.data response_stage4
	$(conv) response_stage4 nic_bend1_004_syn.fits.data response_stage5
	$(conv) response_stage5 nic_pupil_004_syn.fits.data response_stage6

	$(conv) response_stage6 nic1_dqe_008_syn.fits.data response_nic1.res
	$(conv) response_stage6 nic2_dqe_008_syn.fits.data response_nic2.res
	$(conv) response_stage6 nic3_dqe_009_syn.fits.data response_nic3.res



final:
	$(call convolve_nic1,"res")
	$(call convolve_nic2,"res")
	$(call convolve_nic3,"res")

rename:
	rename _012_syn.fits.data.res .flt *
	rename _005_syn.fits.data.res .flt *
	rename nic HST_NICMOS_NIC *.flt
	rename _f _F *.flt
	rename _g _G *.flt
	rename m. M. *.flt
	rename n. N. *.flt
	rename w. W. *.flt
	rename .flt .dat *


calib:
	ls -1 HST*.dat | awk '{ print $$1,$$1,0,0,0; }' >filter.raw
	/work/galev/cocos/calibcosmo /work/galev/filter/calib/vegaabst.dat filter.raw  >filter.dat
	emacs filter.dat

