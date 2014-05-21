function model=stammIdToModel(modelId)
% STAMMIDTOMODEL Convert model ID to model name string.

switch modelId
  case 1
    model='stammIps2State';
  case 2
    model='stammIps3State';
  case 3
    model='stammIps3StateFork';
  case 4
    model='stammIps4StateDark';
  case 5
    model='stammIps4StateFull';
  case 6
    model='stammIps4StateLinear';
  case 7
    model='stammIps4StateFork';
  case 8
    model='stammIps5StateLinear';
  case 9
    model='stammIps4StateCycle';
  case 10
    model='stammIps2StateFwd';
  case 11
    model='stammIps4StateFwd';
  case 12
    model='stammIps5StateFwd';
  case 13
    model='stammIps3StateFwd';
  case 14
    model='stammIps3StateForkFwd';
  case 15
    model='stammIps4StateForkFwd';
  case 16
    model='stammIps4StateDarkFwd';
  case 17
    model='stammIps4StateFwdAna';
  otherwise
    error('Unknown model id');
end
