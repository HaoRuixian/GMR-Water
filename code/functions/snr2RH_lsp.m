function [refl_h, id, psd, pks] = snr2RH_lsp(sinelv, snr, wave_length, hell, hgtlim)
p = polyfit(sinelv, snr, 3);
dsnr = snr - polyval(p, sinelv);

maxf1 = numel(sinelv) / (2*(max(sinelv)-min(sinelv)));
prec1 = 0.001;
ovs = round(wave_length/(2*prec1*(max(sinelv)-min(sinelv))));
fi = 1:1:maxf1;
[psd,f,~,~,~,~] = fLSPw(sinelv,dsnr,fi,0.05,ovs);
psd = cell2mat(psd);
f = cell2mat(f);
refl_h = f.*0.5*wave_length;
surface_h = hell - refl_h;

valid_indx = surface_h>hgtlim(1) & surface_h<hgtlim(2);
refl_h = refl_h(valid_indx);
psd = psd(valid_indx);
[~,id] = max(psd(:));

try
    pks = findpeaks(psd);
catch
    pks = nan;
    return
end
pks = sort(pks);
end