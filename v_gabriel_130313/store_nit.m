function cum_nit_sink = store_nit(ptsnorg, ptscorg, ptup, cupt, pthyd, chyd, ptnit, cnit, ptvol, cvol, ptden, cden,...
 dnorg, dcorg, cberr_norg, cberr_corg,tflnorg,tflcorg,ptmin,cmin, cum_nit_sink)

%Store nitrogen sink terms

cum_nit_sink(end+1,:) = [ptsnorg, ptscorg, ptup', ptmin', pthyd,ptnit,ptvol,ptden, tflnorg, tflcorg, cberr_norg, cberr_corg];
