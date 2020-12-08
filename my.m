data = zeros(8,4);
for i=1:8
    for k=1:4
        ch = ColorHistogramLAB(ct.get_filename_lab(i,k));
        data(i,k) = ch.n_present/ch.n_pixel;
    end
end
