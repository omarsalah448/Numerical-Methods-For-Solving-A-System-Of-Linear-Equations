# VIP, if any of the libraries used here not installed for you
# then open cmd and type "pip install <library name>"
import re
import math
import sympy as sym
from tkinter import *
from utils import *
from tkinter.filedialog import askopenfile

DEFAULT_BACKGROUND = "#424949"
DEFAULT_FONT_COLOR = "#9F1723"
FIELD_BACKGROUND = "#85C1E9"
methodDicNumToWord = {0:"gaussian-elimination",  1:"lu-decomposition", 2:"gaussian-jordan", 3:"gauss-seidel", 4:"jacobi", 5:"all"}
methodDicWordToNum = {"gaussian-elimination":0,  "lu-decomposition":1, "gaussian-jordan":2, "gauss-seidel":3, "jacobi":4, "all":5}
result = 0
methodChoice = 0
b = []
x_initial = []
precision = 9
es = 0.00001
iter_max = 50
eqn = ""
eqns = []
is_file = False

def methodCallback(*args):
    global methodChoice
    for idx,val in methodDicNumToWord.items():
        if val == menu.get():
            methodChoice = idx
    if methodChoice in (0, 1, 2):
        # hide the grid
        x_initalLabel.grid_remove()
        x_initalEntry.grid_remove()
        esLabel.grid_remove()
        esEntry.grid_remove()
        iter_maxLabel.grid_remove()
        iter_maxEntry.grid_remove()
    else:
        # show it again
        x_initalLabel.grid()
        x_initalEntry.grid()
        esLabel.grid()
        esEntry.grid()
        iter_maxLabel.grid()
        iter_maxEntry.grid()

def open_file():
    global is_file, no_eqns, methodChoice, eqns, b, x_initial
    file = askopenfile(mode='r', filetypes=[('Text Files', '*.txt')])
    if file:
        is_file = True
        no_eqns = int(file.readline())
        # remove spaces and new lines then get number from dictionary
        methodChoice = methodDicWordToNum[' '.join(file.readline().split()).lower()]
        eqns = [line.rstrip() for line in file]
        if methodChoice in (3, 4, 5):
            # convert into a list of ints
            x_initial = list(map(int, eqns.pop().split()))
        b = extractAllResults(eqns)
        file.close()

def displayResult():
    global result, b, x_initial, precision, iterations, es, iter_max, eqn, eqns, variables
    if(not is_file):
        try:
            x_initial = list(map(float, x_initalEntry.get().split()))
        except:
            print("ERORR")
        try:
            precision = int(precisionEntry.get())
        except:
            precision = 9
        try:
            es = float(esEntry.get())
        except:
            es = 0.00001
        try:
            iter_max = int(iter_maxEntry.get())
        except:
            iter_max = 50
        try:
            eqn = eqnEntry.get()
        except:
            print("ERORR")
    no_eqns = len(eqns)
    singleStepArray = []
    start_time = time.time()
    variables = extractAllVariables(eqns)
    if methodChoice == 0:
        singleStepArray = gaussian_elimination( eqns )
    elif methodChoice == 1: 
        singleStepArray = LUComp( eqns )
    elif methodChoice == 2: 
        singleStepArray = gauss_jordon( eqns )
    elif methodChoice == 3: 
        singleStepArray = gauss_seidel( eqns, es, iter_max, initial_values=x_initial )
        iterations = [i for i in range(len(singleStepArray))]
        for i in range(len(singleStepArray[0])):
            plt.figure()
            plt.title("Root "+variables[i])
            plot(iterations, column(singleStepArray, i) ,"Iterations", variables[i])
        plt.ylabel("Coefficients");
    elif methodChoice == 4: 
        singleStepArray = jacobi( eqns, es, iter_max, initial_values=x_initial )
        iterations = [i for i in range(len(singleStepArray))]
        for i in range(len(singleStepArray[0])):
            plt.figure()
            plt.title("Root "+variables[i])
            plot(iterations, column(singleStepArray, i) ,"Iterations", variables[i])
        plt.ylabel("Coefficients");
    elif methodChoice == 5: 
        singleStepArray = gaussian_elimination( eqns )
        end_time = time.time()
        outputToFile("w", methodDicNumToWord[0], singleStepArray, no_eqns, precision, (end_time-start_time))
      
        singleStepArray = LUComp( eqns )
        end_time = time.time()
        outputToFile("a", methodDicNumToWord[1], singleStepArray, no_eqns, precision, (end_time-start_time))

        singleStepArray = gauss_jordon( eqns )
        end_time = time.time()
        outputToFile("a", methodDicNumToWord[2], singleStepArray, no_eqns, precision, (end_time-start_time))
        
        singleStepArray = gauss_seidel( eqns, es, iter_max, initial_values=x_initial )
        end_time = time.time()
        outputToFile("a", methodDicNumToWord[3], singleStepArray, no_eqns, precision, (end_time-start_time))

        singleStepArray2 = jacobi( eqns, es, iter_max, initial_values=x_initial )
        end_time = time.time()
        outputToFile("a", methodDicNumToWord[4], singleStepArray, no_eqns, precision, (end_time-start_time))
        iter_seidel = [i for i in range(len(singleStepArray))]
        iter_jacobi = [i for i in range(len(singleStepArray2))]
        for i in range(len(singleStepArray[0])):
            plt.figure()
            plt.title("Root  "+variables[i])
            plot(iter_seidel, column(singleStepArray, i) ,"Iterations", variables[i])
            plot(iter_jacobi, column(singleStepArray2, i) ,"Iterations", variables[i])
            plt.legend(["Gauss-Seidel", "Jacobi"])
        plt.ylabel("Coefficients");
    else:
        result = ""   
    result = "Results outputed to file\n\n"
    answerLabel.config(text=result)
    if methodChoice != 5: 
        end_time = time.time()
        outputToFile("w", methodDicNumToWord[methodChoice], singleStepArray, no_eqns, precision, (end_time-start_time))
        result += "\titer\t"
        for i in range(no_eqns):
            result += variables[i] + '\t'
        result += '\n'
        i = -1
        for tup in singleStepArray:
            i += 1
            result += '\n' + '\t' + str(i) + '\t'
            for elem in tup:
                result += str(elem)[:7] +'\t'
        answerLabel.config(text=result)
def set_text(ch):
    eqnEntry.insert(END,ch)
    return
def del_text():
    eqnEntry.delete(len(eqnEntry.get())-1, END)
    return

def addEqn():
    result = ''
    eqns.append(eqnEntry.get())
    eqnEntry.delete(0, END)
    for eqn in eqns:
        result = result + eqn +'\n'
    matrixLabel.config(text=result)

window = Tk()
window.geometry("1000x500")
# make window responsive
window.columnconfigure(0, weight=1, minsize=75)
window.rowconfigure(0, weight=1, minsize=50)
# create a frame
frame = Frame(width=1000, height=500, background=DEFAULT_BACKGROUND)
frame.pack(anchor=W, fill='both', expand=True)
# title label
titleLabel = Label(master=frame,
                        text="Root Finder",
                        foreground=DEFAULT_FONT_COLOR,
                        background=DEFAULT_BACKGROUND,
                        font=("Times",26,"bold italic"))
titleLabel.grid(row=0, column=11, sticky='n')
# method label
methodLabel = Label(master=frame,
                        text="Method",
                        foreground=DEFAULT_FONT_COLOR,
                        background=DEFAULT_BACKGROUND,
                        font=("Helvetica",26,"bold"))
methodLabel.grid(row=4, column=0, sticky='w')
# method drop down menu
menu = StringVar(frame)
menu.set("Choose A Method")
menu.trace_add('write',methodCallback)
methodMenu= OptionMenu(frame, menu,*methodDicNumToWord.values())
methodMenu.config(background=FIELD_BACKGROUND)
methodMenu.grid(row=4, column=2, sticky='w')
# equation label
eqnLabel = Label(master=frame,
                        text="Equation",
                        foreground=DEFAULT_FONT_COLOR,
                        background=DEFAULT_BACKGROUND,
                        font=("Helvetica",26,"bold"))
eqnLabel.grid(row=6, column=0, sticky='w')
# equation entry
eqnEntry = Entry(master=frame, background=FIELD_BACKGROUND)
eqnEntry.grid(row=6, column=2, sticky='w')
# equation entry ok button
eqnEntryButton = Button(master=frame, text="OK", background=FIELD_BACKGROUND, width=3, command=lambda:addEqn())
eqnEntryButton.grid(row=6, column=3, sticky='w')
# create buttons
buttonOne = Button(master=frame, text="1", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("1"))
buttonOne.grid(row=8, column=0, sticky='w')
buttonTwo = Button(master=frame, text="2", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("2"))
buttonTwo.grid(row=8, column=0)
buttonThree = Button(master=frame, text="3", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("3"))
buttonThree.grid(row=8, column=0, sticky='e')
buttonFour = Button(master=frame, text="4", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("4"))
buttonFour.grid(row=9, column=0, sticky='w')
buttonFive = Button(master=frame, text="5", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("5"))
buttonFive.grid(row=9, column=0)
buttonSix = Button(master=frame, text="6", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("6"))
buttonSix.grid(row=9, column=0, sticky='e')
buttonSeven = Button(master=frame, text="7", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("7"))
buttonSeven.grid(row=10, column=0, sticky='w')
buttonEight = Button(master=frame, text="8", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("8"))
buttonEight.grid(row=10, column=0)
buttonNine = Button(master=frame, text="9", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("9"))
buttonNine.grid(row=10, column=0, sticky='e')
buttonZero = Button(master=frame, text="0", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("0"))
buttonZero.grid(row=11, column=0, sticky='w')
buttonPlus = Button(master=frame, text="+", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("+"))
buttonPlus.grid(row=11, column=0, sticky='e')
buttonMinus = Button(master=frame, text="-", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("-"))
buttonMinus.grid(row=11, column=0)
buttonMul = Button(master=frame, text="*", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("*"))
buttonMul.grid(row=12, column=0, sticky='w')
buttonDiv = Button(master=frame, text="/", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("/"))
buttonDiv.grid(row=12, column=0)
buttonPoint = Button(master=frame, text=".", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("."))
buttonPoint .grid(row=12, column=0, sticky='e')
buttonDel = Button(master=frame, text="Del", background=FIELD_BACKGROUND, width=3, command=del_text)
buttonDel.grid(row=13, column=0, sticky='w')
buttonLeftBracket = Button(master=frame, text="(", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("("))
buttonLeftBracket.grid(row=13, column=0)
buttonRightBracket = Button(master=frame, text=")", background=FIELD_BACKGROUND, width=3, command=lambda:set_text(")"))
buttonRightBracket.grid(row=13, column=0, sticky='e')
buttonVar = Button(master=frame, text="X", background=FIELD_BACKGROUND, width=3, command=lambda:set_text("x"))
buttonVar.grid(row=14, column=0)

# get parameters
parametersLabel = Label(master=frame,
                        text="Parameters",
                        foreground=DEFAULT_FONT_COLOR,
                        background=DEFAULT_BACKGROUND,
                        font=("Helvetica",26,"bold"))
parametersLabel.grid(row=4, column=12, sticky='w')

# get x_inital
x_initalLabel = Label(master=frame,
                text="X Init",
                foreground=DEFAULT_FONT_COLOR,
                background=DEFAULT_BACKGROUND,
                font=("Helvetica",16,"bold"))
x_initalLabel.grid(row=6, column=12, sticky='w')
x_initalEntry = Entry(master=frame, background=FIELD_BACKGROUND)
x_initalEntry.grid(row=6, column=12, sticky='e')

# get precision
precisionLabel = Label(master=frame,
                text="Prec.",
                foreground=DEFAULT_FONT_COLOR,
                background=DEFAULT_BACKGROUND,
                font=("Helvetica",16,"bold"))
precisionLabel.grid(row=7, column=12, sticky='w')
precisionEntry = Entry(master=frame, background=FIELD_BACKGROUND)
precisionEntry.grid(row=7, column=12, sticky='e')

# get es
esLabel = Label(master=frame,
                text="es",
                foreground=DEFAULT_FONT_COLOR,
                background=DEFAULT_BACKGROUND,
                font=("Helvetica",16,"bold"))
esLabel.grid(row=8, column=12, sticky='w')
esEntry = Entry(master=frame, background=FIELD_BACKGROUND)
esEntry.grid(row=8, column=12, sticky='e')
# get iter_max
iter_maxLabel = Label(master=frame,
                text="Iters",
                foreground=DEFAULT_FONT_COLOR,
                background=DEFAULT_BACKGROUND,
                font=("Helvetica",16,"bold"))
iter_maxLabel.grid(row=9, column=12, sticky='w')
iter_maxEntry = Entry(master=frame, background=FIELD_BACKGROUND)
iter_maxEntry.grid(row=9, column=12, sticky='e')
# get file
fileLabel = Label(master=frame,
                text="Choose A File",
                foreground=DEFAULT_FONT_COLOR,
                background=DEFAULT_BACKGROUND,
                font=("Helvetica",16,"bold"))
fileLabel.grid(row=15, column=0, sticky='w')
fileButton = Button(master=frame, text="Browse", background=FIELD_BACKGROUND ,command=open_file)
fileButton.grid(row=16, column=0)

# solve equations button
SEButton = Button(master=frame, text="SOLVE", background=FIELD_BACKGROUND ,command=displayResult)
SEButton.grid(row=11, column=12)
# answer label
answerTextLabel = Label(master=frame,
                        text="ANSWER:",
                        foreground=DEFAULT_FONT_COLOR,
                        background=DEFAULT_BACKGROUND,
                        font=("Times",26,"bold italic"))
answerTextLabel.grid(row=14, column=11, sticky='n')
answerLabel = Label(master=frame,
                        text="",
                        foreground="cyan",
                        background=DEFAULT_BACKGROUND,
                        font=("Helvetica",12,"bold"))
answerLabel.grid(row=14, column=12, sticky='w')
matrixLabel = Label(master=frame,
                        text="",
                        foreground="cyan",
                        background=DEFAULT_BACKGROUND,
                        font=("Helvetica",12,"bold"))
matrixLabel.grid(row=7, column=2)
window.mainloop()
