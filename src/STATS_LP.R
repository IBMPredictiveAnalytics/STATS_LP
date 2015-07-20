#Licensed Materials - Property of IBM
#IBM SPSS Products: Statistics General
#(c) Copyright IBM Corp. 2014
#US Government Users Restricted Rights - Use, duplication or disclosure 
#restricted by GSA ADP Schedule Contract with IBM Corp.

# Author: JKP, IBM SPSS
# Version = 1.0.1

# history
# 17-Jun-2013 - original version
# 05-Jul-2013 - rewrite to base on package lpSolveAPI instead of linprog


helptext="STATS LP OBJECTIVE = list of objective function values
CBOUNDS = bounds variable name CBOUNDSLBLS = bounds labels
CONSTRAINTS = list of constraint variables
CONSTRAINTDIR = LE, EQ or GE CONSTRAINTDIRVAR = variable name
INTEGERVARS = variable list BINARYVARS = variable list

/OBJBOUNDS objective variable lower and upper bounds

/OPTIONS MAXITER = integer TOLERANCE = positive number
ASZERO = positive number

/OUTPUT CONSTRAINTS = YES* or NO

/HELP

Example:
data list list/x1 x2 x3 b(4F8.2) blbls(a10).
begin data.
.7 .35 0 40 'first'
1.5 1 3 90  'second'
50 12.5 20 2500 'third'
end data.
dataset name lp.
stats lp objective=-1800 -600 -600 cbounds=b cboundslbls=blbls
constraints=x1 x2 x3 constraintdir=le.

OBJECTIVE (required) specifies the coefficients for the objective
function.  They must be listed in the same order as the 
names listed for constraint variables.  There must be one value
for each variable listed under CONSTRAINTS.

CONSTRAINTS (required) lists the variables holding the constraint
values.

CBOUNDS (required) specifies the variable holding the bounds for
each constraint

CONSTRAINTDIR or CONSTRAINTDIRVAR specify the type of
the constraints.  CONSTRAINTDIR can be LE, EQ, or GE indicating
that all constraints are <=, =, or >= the bound.  Alternatively,
CONSTRAINTDIRVAR specifies a variable whose value for each
case (constraint) is '<=', '=', or '>=', allows the directions to differ.  
(There is no difference between < and <=.)  
Either CONSTRAINTDIR or CONSTRAINTDIRVAR must be specified
but not both.

INTEGERVARS and BINARYVARS can specify objective variables
whose values must be integer or 0,1, respectively.  The variables
must also appear in the CONSTRAINTS list.  A variable should
not be included in both lists.  

/OBJBOUNDS
This subcommand is used to specify lower and upper bounds
on the objective variables.  Each entry should have the
form
variable name lower bound upperbound, e.g., X 5 10.
A value written as a period means that that bound
is not specified.  For example,
X . 10
means an upper bound of 10 for X but no lower bound is
specified.
The name-value-value sequence can be repeated as many time as necessary.
Specifying a variable as integer and bounds of 0 and 1 is equivalent
to specifying it as binary.  Changing the bounds on a binary variable
is equivalent to specifying it as integer.

/OPTIONS

MAXITER (optional) specifies the maximum number of iterations and
defaults to 1000.

TOLERANCE (optional) specifies how much a constraint can be
violated before declaring failure.  It defaults to 10**-6.

ASZERO specifies a value below which values are set to
exactly zero.  It defaults to 10**-9.

/OUTPUT

CONSTRAINTS specifies whether to print details about the constraints
at the solution.

STATS LP /HELP prints this help and does nothing else.
"

lpext <- function(objective, objectivedir="min",cbounds, 
	constraints, constraintdir=NULL, 
    constraintdirvar=NULL, cboundslbls=NULL, integervars=NULL, binaryvars=NULL,
    allinteger=FALSE, allbinary=FALSE,
    objdirection="minimize", maxiter=1000, tolerance=1e-6, 
    aszero=1e-9, uselpsolve=FALSE,
    constraintsds=TRUE, allvars=FALSE, objbounds=NULL) {

    setuplocalization("STATS_LP")	

    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Linear Programming")
    warningsprocname = gtxt("Linear Programming: Warnings")
    omsid="STATSLP"
    warns = Warn(procname=warningsprocname,omsid=omsid)

    tryCatch(library(lpSolveAPI), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.","lpSolveAPI"),
            dostop=TRUE)
    }
    )
    if (!xor(is.null(constraintdir), is.null(constraintdirvar))) {
        warns$warn(gtxt("Either a constraint direction or a constraint direction variable must be specified but not both"), dostop=TRUE)
    }
	if (allinteger && allbinary) {
		warns$warn(gtxt("Cannot specify both ALLINTEGERS=YES and ALLBINARY=YES"),
			dostop=TRUE)
	}
    if (allinteger && !is.null("integervars")) {
        warns$warn(gtxt("Cannot specify both a list of integer variables and ALLINTEGER=YES"),
            dostop=TRUE)
    }
    if (allinteger) {
        integervars = constraints
    }
	if (allbinary && !is.null("binaryvars")) {
	warns$warn(gtxt("Cannot specify both a list of binary variables and ALLBINARY=YES"),
		dostop=TRUE)
    }
    if (allbinary) {
        binaryvars = constraints
    }

	if (!is.null(constraintdir)) {
		if (constraintdir == "le") {
			constraintdir = "<="
		} else if (constraintdir == "ge") {
			constraintdir = ">="
		} else if (constraintdir == "eq") { 
			constraintdir = "="
		}
	}
    
    n = length(objective)  # size of objective function
    nx = length(constraints) # number of constraint variables

    if (n != nx) {
        warns$warn(gtxt("The number of objective coefficients differs from the number of constraint variables"), 
            dostop=TRUE)
    }

    alldata = c(constraints, cbounds, constraintdirvar, cboundslbls)
    dta = spssdata.GetDataFromSPSS(alldata, missingValueToNA=TRUE, factorMode="none")
    nconstraints = nrow(dta)
    # blanks must be trimmed from const.dir if used, and const.dir=NULL is not allowed
    bvec=dta[[nx+1]]


    if (!is.null(constraintdirvar)) {
        dta[nx+2] = gsub(" ", "", dta[[nx+2]])
		# lpSolve treats an invalid value as <=.  Check all values for legitimacy
		if (any(is.na(match(unlist(dta[nx+2]), c("<=", "=", ">="))))) {
			warns$warn(gtxt("An invalid constraint direction variable value was found.  Values must be <=, =, or >=."), dostop=TRUE)
		}
    }

    # set up the lp model
    lpmodel = make.lp(nconstraints, nx, verbose="neutral")
	lp.control(lpmodel, sense=objectivedir)
    for (i in 1:nx) {
        set.column(lpmodel, i, dta[,i])
    }
    set.objfn(lpmodel, objective)
	if (!is.null(integervars)) {
		intvars = match(integervars, constraints)
		if (any(is.na(intvars))) {
			warns$warn(gtxt("A variable specified as type integer does not appear in the list of constraint variables"), dostop=TRUE)
		}
		set.type(lpmodel, intvars, "integer")
	}	
	if (!is.null(binaryvars)) {
		binvars = match(binaryvars, constraints)
		if (any(is.na(binvars))) {
			warns$warn(gtxt("A variable specified as type binary does not appear in the list of constraint variables"), dostop=TRUE)
		}
		set.type(lpmodel, binvars, "binary")
	}
    if (!is.null(constraintdir)) {
        set.constr.type(lpmodel, rep(constraintdir, nconstraints))
    } else {
        set.constr.type(lpmodel, dta[[nx+2]])
    }
    set.rhs(lpmodel, bvec)
    if (!is.null(cboundslbls)) {
        dimnames(lpmodel)[[1]] = dta[[length(dta)]]
    }
    dimnames(lpmodel)[[2]] = constraints
    
    setobjbounds(lpmodel, objbounds, warns)
    res = tryCatch(solve(lpmodel),
        error=function(e) {warns$warn(e$message, dostop=TRUE)
            }
        )


    # print results
    # 
    displayresults(lpmodel, objective, res, constraints, cbounds, constraintdir, 
        constraintdirvar,tolerance, aszero,
        constraintsds, allvars, warns)
	# calling delete.lp sometimes causes an error complaining about the model
	# being null, because the model has already had its finalizer called.
    ###delete.lp(lpmodel)    # free external resources
    # clean up workspace
    res <- tryCatch(rm(list=ls()), warning = function(e) {return(NULL)})
}


rounder <- function(s) {
    if (!is.null(s)) {
        s = round(s, 5)
    } else {
        s = "."
    }
    return(s)
}

tamenumber <- function(x) {
	# return x or a string representation in scientific notation
	if (abs(x) < 1e8) {
		return(rounder(x))
	} else {
		return(sprintf("%.4e", x))
	}
}

tamedf = function(df) {
    # modify df to force extreme values into scientific notation strings
	# there must be an easier way to do this!
	return(data.frame(matrix(lapply(unlist(df), tamenumber), nrow(df), ncol(df)),
		row.names=row.names(df)))
}

setobjbounds = function(lpmodel, objbounds, warns) {
	# set lower and upper bounds for objective variables, if any

	if (!is.null(objbounds)) {
		objlen = length(objbounds)
		if (objlen %% 3 != 0) {  # items have the form varname lowerbound upperbound
			###delete.lp(lpmodel)  leave this to the finalizer
			warns$warn(gtxt("The objbounds statement is missing one or more elements"),
				dostop=TRUE)
		}
		for (i in seq(1, objlen, 3)) {
			varindex = match(objbounds[i], colnames(lpmodel))
			if (is.na(varindex)) {
				delete.lp(lpmodel)		
				warns$warn(gtxtf("A variable was specified on OBJBOUNDS that is  not in the model: %s",
					objbounds[i]), dostop=TRUE)
			}
			if (objbounds[i+1] != ".") {  # "." means no bound
				set.bounds(lpmodel, lower=objbounds[i+1], columns=varindex)
			}
			if (objbounds[i+2] != ".") {
				set.bounds(lpmodel, upper=objbounds[i+2], columns=varindex)
			}
		}
	}
}

displayresults = function(lp, objective, res, constraints, cbounds, constraintdir, 
    constraintdirvar,tolerance, aszero, 
    constraintsds, allvars, warns) {

    StartProcedure(gtxt("Linear Programming"), "STATSLP")
    retco = list(       # return codes
		gtxt("optimal solution found"),				#0
		gtxt("the model is sub-optimal"),           #1 
		gtxt("the model is infeasible"),            #2 
		gtxt("the model is unbounded"),             #3 
		gtxt("the model is degenerate"),            #4 
		gtxt("numerical failure encountered"),      #5 
		gtxt("process aborted"),                    #6 
		gtxt("timeout"),                            #7 
		"unused",                                         
		gtxt("the model was solved by presolve"),   #9 
		gtxt(	"the branch and bound routine failed"), #10
		gtxt(	"the branch and bound was stopped because of a break-at-first or break-at-value"), #11
		gtxt(	"a feasible branch and bound solution was found"),		#12
		gtxt(	"no feasible branch and bound solution w13:as found")  #13
	)
    
    # summary results
    scaption = gtxt("Computations done by R package lpSolve")
    lbls = c(
        gtxt("Objective Function Optimal Value"), 
        gtxt("Solution Status"), 
		gtxt("Objective Direction"),
        gtxt("Constraint Variables"),
        gtxt("Bounds"), 
        gtxt("Constraint Direction"), 
        gtxt("Total Iterations"), 
        gtxt("Iterations Used for Relaxed Solution"), 
        gtxt("Iterations Used for Branch and Bound"), 
        gtxt("Tolerance"), 
        gtxt("Treat As Zero")
    )
	niter = get.total.iter(lp)
	if (length(niter) == 1) {
		totaliter = niter
		hasintegers = FALSE
	} else {
		totaliter = niter[1]+niter[2]
		hasintegers = TRUE
	}
    vals = c(
        rounder(get.objective(lp)),
        retco[res+1],
		ifelse(lp.control(lp)$sense == "minimize", gtxt("minimize"), gtxt("maximize")),
        paste(constraints, collapse = " "),
        cbounds,
        ifelse(is.null(constraintdir), constraintdirvar, constraintdir),
        totaliter,
        ifelse(hasintegers, niter[1], gtxt("--NA--")),
        ifelse(hasintegers, niter[2], gtxt("--NA--")),
        lp.control(lp)$epsilon[3],
        lp.control(lp)$epsilon[1]
    )

    # settings and result summary
    spsspivottable.Display(data.frame(cbind(vals), row.names=lbls), title = gtxt("Settings and Results Summary"),
        collabels=c(gtxt("Summary")), templateName="LPSUMMARY", outline=gtxt("Summary"),
        caption = scaption)
      
    # objective function
	# f maps unbounded elements to "None"
	f = function(item) {if (item == Inf) {return(gtxt("None"))} else {return(item)}}
	bnds = get.bounds(lp)
	bl = lapply(bnds$lower, f)
	bu = lapply(bnds$upper, f)

    spsspivottable.Display(data.frame(objective, get.type(lp), cbind(bl, bu), row.names=constraints),
        title=gtxt("Objective Function Coefficients and Bounds"), collabels=c(gtxt("Coefficients"),
			gtxt("Variable Type"), gtxt("Lower Bound"), gtxt("Upper Bound")),
        templateName="LPOBJECTIVE", outline=gtxt("Objective Function"))

    if (res == 0) {
        
        # optimal solution variables
        spsspivottable.Display(data.frame(get.variables(lp), row.names=constraints), title=gtxt("Variables at Optimum"),
            collabels=c(gtxt("Values")), templateName="LPBASIS", outline=gtxt("Optimum"),
		)

		#objective function sensitivity
		so = data.frame(get.sensitivity.obj(lp), row.names=constraints)
		so = tamedf(so)
		names(so) = c(gtxt("From"), gtxt("To"))
		spsspivottable.Display(so, title=gtxt("Objective Function Sensitivity"),
			templateName="LPOBJSENS",
			outline=gtxt("Objective Function Sensitivity"),
			caption=gtxt("Rows show the range within which that objective function coefficient can vary without changing the solution.  However, they do not account for any integer constraints. ")
		)
        # constraint usage
		bound = get.constr.value(lp, side="rhs")
		lbound = length(bound)
		actual = get.constraints(lp)
		dual = get.dual.solution(lp)
		shadow = -dual[(1:lbound)+1]
		shadowbnds = get.sensitivity.rhs(lp)
		consdf = data.frame(actual, get.constr.type(lp), bound, bound-actual, shadow,
			shadowbnds$dualsfrom[1:lbound], shadowbnds$dualstill[1:lbound],
			row.names=dimnames(lp)[[1]])
        if (constraintsds) {
			consdf = tamedf(consdf)
			spsspivottable.Display(
				consdf,
				title=gtxt("Constraints at Optimum"),
				collabels = c(gtxt("Actual"), gtxt("Direction"), gtxt("Bound"), 
				gtxt("Difference"), gtxt("Shadow Price"), 
				gtxt("Lower Bound"), gtxt("Upper Bound")),
				templateName="LPCONSTRAINTS",
				caption=gtxt("Lower and upper bounds show the range within which the bound can vary without changing the solution.  However, they do not account for any integer constraints. "),
				outline=gtxt("Constraint Results")
			)
        }
    }
    warns$display(inproc=TRUE)
}

# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}
gtxt <- function(...) {
    return(gettext(...,domain="STATS_LP"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_LP"))
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

    if (is.null(msg) || dostop) {
        lcl$display(inproc)  # display messages and end procedure state
        if (dostop) {
            stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
        }
    }
}

    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

    if (lcl$msgnum == 0) {   # nothing to display
        if (inproc) {
            spsspkg.EndProcedure()
        }
    } else {
        if (!inproc) {
            procok =tryCatch({
                StartProcedure(lcl$procname, lcl$omsid)
                TRUE
                },
                error = function(e) {
                    FALSE
                }
            )
        }
        if (procok) {  # build and display a Warnings table if we can
            table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
            rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

    for (i in 1:lcl$msgnum) {
        rowcategory = spss.CellText.String(as.character(i))
        BasePivotTable.SetCategories(table,rowdim,rowcategory)
        BasePivotTable.SetCellValue(table,rowcategory, 
            spss.CellText.String(lcl$msglist[[i]]))
    }
    spsspkg.EndProcedure()   # implies display
} else { # can't produce a table
    for (i in 1:lcl$msgnum) {
        print(lcl$msglist[[i]])
    }
}
}
}
return(lcl)
}

Run<-function(args){
    cmdname = args[[1]]
    args <- args[[2]]
    oobj<-spsspkg.Syntax(templ=list(
        spsspkg.Template("OBJECTIVE", subc="",  ktype="float", var="objective", islist=TRUE),
		spsspkg.Template("OBJECTIVEDIR", subc="", ktype="str", var="objectivedir",
			vallist=list("min", "max")),
        spsspkg.Template("CBOUNDS", subc="",  ktype="existingvarlist", var="cbounds"),
        spsspkg.Template("CBOUNDSLBLS", subc="", ktype="existingvarlist", var="cboundslbls"),
        spsspkg.Template("CONSTRAINTS", subc="",  ktype="existingvarlist", 
            var="constraints", islist=TRUE),
        spsspkg.Template("CONSTRAINTDIR", subc="",  ktype="str", 
            var="constraintdir", vallist=list("le", "ge", "eq")),
        spsspkg.Template("CONSTRAINTDIRVAR", subc="", ktype="existingvarlist",
            var="constraintdirvar"),
        spsspkg.Template("INTEGERVARS", subc="", ktype="existingvarlist", 
            var="integervars", islist=TRUE),
		spsspkg.Template("BINARYVARS", subc="", ktype="existingvarlist", 
            var="binaryvars", islist=TRUE),
        spsspkg.Template("ALLINTEGER", subc="", ktype="bool", var="allinteger"),
		spsspkg.Template("ALLBINARY", subc="", ktype="bool", var="allbinary"),

        spsspkg.Template("MAXITER", subc="OPTIONS", ktype="int", 
            var="maxiter", vallist = list(1)), 
        spsspkg.Template("TOLERANCE", subc="OPTIONS", ktype="float", 
            var="tolerance", vallist=list(1e-20)),
        spsspkg.Template("ASZERO", subc="OPTIONS", ktype="float", 
        var="aszero", vallist=list(0)),		
        
        spsspkg.Template("CONSTRAINTS", subc="OUTPUT", ktype="bool", var="constraintsds"),
        spsspkg.Template("ALLVARS", subc="OUTPUT", ktype="bool", var="allvars"),
		
		spsspkg.Template("", subc="OBJBOUNDS", ktype="literal", var="objbounds", islist=TRUE)
        ))        
    if ("HELP" %in% attr(args,"names")) {
        #writeLines(helptext)
        helper(cmdname)
    } else {
        res <- spsspkg.processcmd(oobj,args,"lpext")
    }
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}