VERSION 5.00
Begin VB.Form Form1 
   Caption         =   "VBTemplate"
   ClientHeight    =   3192
   ClientLeft      =   60
   ClientTop       =   348
   ClientWidth     =   4680
   LinkTopic       =   "Form1"
   PaletteMode     =   1  'UseZOrder
   ScaleHeight     =   3192
   ScaleWidth      =   4680
   StartUpPosition =   3  'Windows Default
   Begin VB.TextBox Text1 
      Height          =   855
      Left            =   840
      TabIndex        =   1
      Text            =   "Text1"
      Top             =   600
      Width           =   1215
   End
   Begin VB.CommandButton RunModel 
      Caption         =   "Run "
      Height          =   855
      Left            =   3240
      TabIndex        =   0
      Top             =   720
      Width           =   975
   End
End
Attribute VB_Name = "Form1"
Attribute VB_GlobalNameSpace = False
Attribute VB_Creatable = False
Attribute VB_PredeclaredId = True
Attribute VB_Exposed = False
Option Explicit
'following are global variables

Private Sub form_load()
'this gets called before the click procedure
Dim time As Single
Dim I As Integer
Dim J As Integer
Dim Q As Single
Open "shellboundary.csv" For Input As #1
    Input #1, time
Close #1

Open "boundary.csv" For Output As #1
    Print #1, "SV  Compartment Boundary Value(g/m3)"
 
    For I = 1 To 8
        Q = 2000 ' default
        Print #1, "1,"; I; ",";
        If I = 1 Then
            If time < 3600 And time >= 720 Then
                Q = 2000 + 1.389 * (time - 720)
            End If
            If time < 6480 And time >= 3600 Then
                Q = 6000 - 1.389 * (time - 3600)
            End If
            Print #1, Q
        End If
        If I > 1 And I < 8 Then
            Print #1, -999
        End If
        If I = 8 Then
            Print #1, 2000
        End If
    Next I
Close #1
End
End Sub



