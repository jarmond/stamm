function id=stammModelToId(model)
% STAMMMODELTOID Convert model string to numerical ID

switch model
  case 'stammIps2State'
    id=1;
  case 'stammIps3State'
    id=2;
  case 'stammIps3StateFork'
    id=3;
  case 'stammIps4StateDark'
    id=4;
  case 'stammIps4StateFull'
    id=5;
  case 'stammIps4StateLinear'
    id=6;
  case 'stammIps4StateFork'
    id=7;
  case 'stammIps5StateLinear'
    id=8;
  case 'stammIps4StateCycle'
    id=9;
  case 'stammIps2StateFwd'
    id=10;
  case 'stammIps4StateFwd'
    id=11;
  case 'stammIps5StateFwd'
    id=12;
  case 'stammIps3StateFwd'
    id=13;
  case 'stammIps3StateForkFwd'
    id=14;
  case 'stammIps4StateForkFwd'
    id=15;
  case 'stammIps4StateDarkFwd'
    id=16;
  case 'stammIps4StateFwdAna'
    id=17;
  otherwise
    error('Unknown model');
end
