function tblOut = tblReorg(EEG_table, Spec_table, Epoch_table, EpochSpec_table)

tblOut = table();


% baselines --------------------------------------------------------------

tblOut.BaselineOpenBefore('EEG') = EEG_table.BaselineOpen('before experiment');
tblOut.BaselineOpenAfter('EEG')  = EEG_table.BaselineOpen('after experiment');
tblOut.BaselineClosedBefore('EEG') = EEG_table.BaselineClosed('before experiment');
tblOut.BaselineClosedAfter('EEG')  = EEG_table.BaselineClosed('after experiment');

tblOut.BaselineOpenBefore('Spec') = Spec_table.BaselineOpen('before experiment');
tblOut.BaselineOpenAfter('Spec')  = Spec_table.BaselineOpen('after experiment');
tblOut.BaselineClosedBefore('Spec') = Spec_table.BaselineClosed('before experiment');
tblOut.BaselineClosedAfter('Spec')  = Spec_table.BaselineClosed('after experiment');

tblOut.BaselineOpenBefore('Epoch EEG') = Epoch_table.BaselineOpen('before experiment');
tblOut.BaselineOpenAfter('Epoch EEG')  = Epoch_table.BaselineOpen('after experiment');
tblOut.BaselineClosedBefore('Epoch EEG') = Epoch_table.BaselineClosed('before experiment');
tblOut.BaselineClosedAfter('Epoch EEG')  = Epoch_table.BaselineClosed('after experiment');

tblOut.BaselineOpenBefore('Epoch Spec') = EpochSpec_table.BaselineOpen('before experiment');
tblOut.BaselineOpenAfter('Epoch Spec')  = EpochSpec_table.BaselineOpen('after experiment');
tblOut.BaselineClosedBefore('Epoch Spec') = EpochSpec_table.BaselineClosed('before experiment');
tblOut.BaselineClosedAfter('Epoch Spec')  = EpochSpec_table.BaselineClosed('after experiment');


% stims ------------------------------------------------------------------

tblOut.TempStim('EEG') = UnwrapCat(EEG_table.TempStim('before experiment'), EEG_table.TempStim('after experiment'));
tblOut.PinPrick('EEG') = UnwrapCat(EEG_table.PinPrick('before experiment'), EEG_table.PinPrick('after experiment'));
tblOut.Pressure('EEG') = UnwrapCat(EEG_table.Pressure('before experiment'), EEG_table.Pressure('after experiment'));

tblOut.TempStim('Spec') = UnwrapCat(Spec_table.TempStim('before experiment'), Spec_table.TempStim('after experiment'));
tblOut.PinPrick('Spec') = UnwrapCat(Spec_table.PinPrick('before experiment'), Spec_table.PinPrick('after experiment'));
tblOut.Pressure('Spec') = UnwrapCat(Spec_table.Pressure('before experiment'), Spec_table.Pressure('after experiment'));

tblOut.TempStim('Epoch EEG') = UnwrapCat(Epoch_table.TempStim('before experiment'), Epoch_table.TempStim('after experiment'));
tblOut.PinPrick('Epoch EEG') = UnwrapCat(Epoch_table.PinPrick('before experiment'), Epoch_table.PinPrick('after experiment'));
tblOut.Pressure('Epoch EEG') = UnwrapCat(Epoch_table.Pressure('before experiment'), Epoch_table.Pressure('after experiment'));

tblOut.TempStim('Epoch Spec') = UnwrapCat(EpochSpec_table.TempStim('before experiment'), EpochSpec_table.TempStim('after experiment'));
tblOut.PinPrick('Epoch Spec') = UnwrapCat(EpochSpec_table.PinPrick('before experiment'), EpochSpec_table.PinPrick('after experiment'));
tblOut.Pressure('Epoch Spec') = UnwrapCat(EpochSpec_table.Pressure('before experiment'), EpochSpec_table.Pressure('after experiment'));


% CPM --------------------------------------------------------------------

tblOut.PinPrickCPM('EEG') = EEG_table.PinPrick('CPM');
tblOut.PressureCPM('EEG') = EEG_table.Pressure('CPM');

tblOut.PinPrickCPM('Spec') = Spec_table.PinPrick('CPM');
tblOut.PressureCPM('Spec') = Spec_table.Pressure('CPM');

tblOut.PinPrickCPM('Epoch EEG') = Epoch_table.PinPrick('CPM');
tblOut.PressureCPM('Epoch EEG') = Epoch_table.Pressure('CPM');

tblOut.PinPrickCPM('Epoch Spec') = EpochSpec_table.PinPrick('CPM');
tblOut.PressureCPM('Epoch Spec') = EpochSpec_table.Pressure('CPM');

end