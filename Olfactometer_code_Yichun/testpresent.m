function testpresent(Seq)

Duration=1;
ISI=25;
Pre=5;

ISI=ISI-Pre;

FlipValve_EP('all',1)

for i=1:length(Seq)
    injectOdour_EP(Seq(i))
    pause(Pre)
    tic
    FlipValve_EP('Final',0)
    pause(Duration)
    FlipValve_EP('all',1)
    if i<length(Seq)
        pause(ISI-toc)
    end
end
