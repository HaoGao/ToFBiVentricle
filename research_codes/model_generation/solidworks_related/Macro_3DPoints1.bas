Attribute VB_Name = "Macro_3DPoints1"
' ******************************************************************************
' C:\Users\haogao\AppData\Local\Temp\swx3512\Macro1.swb - macro recorded on 04/05/19 by haogao
' ******************************************************************************
Dim swApp As Object

Dim Part As Object
Dim boolstatus As Boolean
Dim longstatus As Long, longwarnings As Long

Sub main()

Set swApp = _
Application.SldWorks

Set Part = swApp.ActiveDoc
Dim myModelView As Object
Set myModelView = Part.ActiveView
myModelView.FrameState = swWindowState_e.swWindowMaximized




Part.SketchManager.Insert3DSketch True
Dim pointArray As Variant
Dim points() As Double
ReDim points(0 To 44) As Double
points(0) = 0.02272
points(1) = 0.13807
points(2) = 0.86194
points(3) = 0.02466
points(4) = 0.13978
points(5) = 0.85501
points(6) = 0.02461
points(7) = 0.13932
points(8) = 0.85256
points(9) = 0.02409
points(10) = 0.13876
points(11) = 0.85034
points(12) = 0.02336
points(13) = 0.13774
points(14) = 0.84948
points(15) = 0.02317
points(16) = 0.1374
points(17) = 0.84717
points(18) = 0.02211
points(19) = 0.13606
points(20) = 0.84551
points(21) = 0.02119
points(22) = 0.13482
points(23) = 0.84446
points(24) = 0.01912
points(25) = 0.13256
points(26) = 0.84445
points(27) = 0.0173
points(28) = 0.13052
points(29) = 0.84398
points(30) = 0.01555
points(31) = 0.12903
points(32) = 0.84508
points(33) = 0.01312
points(34) = 0.12704
points(35) = 0.84765
points(36) = 0.01187
points(37) = 0.1259
points(38) = 0.84856
points(39) = 0.01093
points(40) = 0.12491
points(41) = 0.85069
points(42) = 0.00997
points(43) = 0.12434
points(44) = 0.85433
pointArray = points
Set skSegment = Part.SketchManager.CreateSpline((pointArray))
Part.ClearSelection2 True
Part.SketchManager.InsertSketch True






End Sub

