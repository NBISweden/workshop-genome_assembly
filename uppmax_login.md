# How to log in to Uppmax 

## Prerequisites:

* A Shell terminal
* The Unix tools `ssh` and `scp`.

If you are on Linux: A terminal should be installed.

If you are on MacOS: A terminal should be installed, however [XQuartz](https://www.xquartz.org/) 
must be installed to run graphical applications.

If you are on Windows: Please download and install [MobaXterm](http://mobaxterm.mobatek.net). 

## 1. Connecting to UPPMAX

The first step of this lab is to open a ssh connection to UPPMAX. To do this you need 
to start a terminal program on your computer. Which program is used depends on the 
system of your computer (MacOSX, Windows or Linux). On MacOSX you should use Terminal 
(icon looks like a black screen), on Windows you should use 
[MobaXterm](http://mobaxterm.mobatek.net). For Linux it depends on which version of 
Linux you are running, but it should be fairly obvious which program to use.

Once you have started your terminal program, type in (change `<username>` to your own username):

```bash
$ ssh -X <username>@rackham.uppmax.uu.se
```

and give your password when prompted. As you type the password, nothing will show on 
screen. No stars, no dots. It is supposed to be that way. Just type the password and 
press enter.

You should now get a welcoming message from Uppmax to show that you have successfully 
logged in.

## 2. Logging into the reserved node

For this course, we have arranged for you to have one half of a node (=10 cores) each. 
To get this reservation you need to use the salloc command like this:

```bash
$ salloc -A <reservation> -t 08:00:00 -p core -n 10 --no-shell --reservation=<reservation_id> &
```

where you should substitute `<reservation>` with `` and `<reservation_id>` with one of these 
alternatives depending on the day.

Day 1: `g2018025_1`

Day 2: `g2018025_2`

Day 3: `g2018025_3`

Now check which node you got (replace *username* with your uppmax user name) like this:

```bash
$ squeue -u username
```

Under Nodelist you will see the name of the node that has been reserved for you. The names 
follow the format mXX, e.g., m34.

Connect to the node you were dealt like this:

```bash
$ ssh -X m34
```

**Note**: there is a uppmax specific tool called jobinfo that supplies the same kind of 
information as squeue that you can use as well ( `$ jobinfo -u username`).

You are now logged in to your reserved node, and there is no need for you to use the 
SLURM queuing system. You can now continue with the specific exercise instructions, 
but **please remember to never use more than 10 cores** as you will be sharing a node 
(i.e., a single computer) with someone else in the course and there are only 20 cores 
available in total per node.

## Helpful Resources

[Linux cheat sheet](files/linux-cheat-sheet.pdf)

[UPPMAX cheat sheet](files/uppmax-cheat-sheet.png)

