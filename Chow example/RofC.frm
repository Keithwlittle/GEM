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
Dim RofC(6) As Single
Dim X(6) As Single
Dim alpha As Single
Dim beta As Single
alpha = 3.49 '* 86400 / 3.28 ' convert from sec/ft to day/m
beta = 0.6
Open "shellRofC.csv" For Input As #1
    Input #1, time
    For I = 1 To 6
        Input #1, X(I)
    Next I
Close #1

Open "RofC.csv" For Output As #1
For I = 1 To 6
    RofC(I) = alpha * X(I) ^ (beta - 1)
    Print #1, RofC(I)
Next I
Close #1
End
End Sub



