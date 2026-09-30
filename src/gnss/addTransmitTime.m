function raw = addTransmitTime(raw, c)
% ADDTRANSMITTIME 由接收时刻与原始伪距计算信号发射时刻 txTime = rxTime - prRaw/c（GPST, s）。
transmit = num2cell([raw.rxTime] - [raw.prRaw]/c);
[raw.txTime] = transmit{:};
end
