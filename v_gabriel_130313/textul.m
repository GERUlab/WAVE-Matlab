function h = textul(txt,margin1,margin2)
    a = axis;
    wdth = a(2)-a(1);
    ht = a(4)-a(3);
    pos = [a(1)+margin1*wdth a(4)-margin2*ht];
    h = text(pos(1),pos(2),txt);
end