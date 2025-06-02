function w = narrow_focus(Hc,Nt)
w = (Hc./abs(Hc))'*sqrt(1/Nt);
end

