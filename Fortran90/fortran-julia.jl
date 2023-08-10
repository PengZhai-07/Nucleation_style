#=
This julia script converts fortran 90 code into julia.
It uses naive regex replacements to do as much as possible,
but the output WILL need further cleanup.
Known conversion problems such as GOTO are commented and marked with FIXME
Most variable declaration lines are entirely deleted, which may or
may not be useful. 
To run from a shell: 
julia fortran-julia.jl filename.f90
Output is written to filename.jl.
=#

using DataStructures

# Regex/substitution pairs for replace(). Order matters here.
replacements = OrderedDict(
  # Lowercase everything not commented
  r"^(?!.*!).*"m => lowercase,
  # Lowercase start of lines with comments
  r"^.*!"m => lowercase,
  # Remove '&' multiline continuations
  r"\s*&\s*" => "",
  # Comments use # not !
  "!" => "#",
  # Powers use ^ not **
  "**" => "^",
  # Only double quotes allowed for strings
  "'" => "\"",
  # DO loop to for loop
  r"do (.*),(.*)" => s"for \1:\2",
  # Spaces around math operators
  r"([\*\+\/=])(?=\S)" => s"\1 ",
  r"(?<=\S)([\*\+\/=])" => s" \1",
  # Spaces around - operators, except after e
  # r"([^e][\-])(\S)" => s"\1 \2",
  r"(?<!\W\de)(\h*\-\h*)" => s" - ",
  # Space after all commas
  r"(,)(\S)" => s"\1 \2",
  # Replace ELSEIF/ELSE IF with elseif 
  r"(\s+)else if" => s"\1elseif",
  # Replace IF followed by ( to if (
  r"(\s+)(elseif|if)\(" => s"\1\2 (",
  # Remove THEN
  r"([)\s])then(\s+)" => s"\1\2",
  # Relace END XXXX with end
  r"(\s+)end\h*.*" => s"\1end",
  # Replace expnent function
  r"(\W)exp\(" => s"\1exp(",
  # Reorganise functions and doc strings. This may be very project specific.
  r"#\^\^+\s*subroutine\s*(\w+)([^)]+\))\s*(.*?)#\^\^\^+"sm => 
      Base.SubstitutionString("\"\"\"\n\\3\"\"\"\nfunction \\1\\2::Void"),
  r"\#\^\^+\s*real function\s*(\w+)([^)]+\))\s*(.*?)\#\^\^\^+"sm => 
      Base.SubstitutionString("\"\"\"\n\\3\"\"\"\nfunction \\1\\2::Float64"),
  # Don't need CALL
  r"(\s*)call(\h+)" => s"\1",
  # Use real math symbols
  "gamma" => "Γ",
  "theta" => "Θ",
  "epsilon" => "ϵ",
  "lambda" => "λ",
  "alpha" => "α",
  # Swap logical symbols
  ".true." => "true",
  ".false." => "false",
  r"\s*\.or\.\s*" => " || ",
  r"\s*\.and\.\s*" => " && ",
  r"\s*\.not\.\s*" => " ! ",
  r"\s*\.eq\.\s*" => " == ",
  r"\s*\.ne\.\s*" => " != ",
  r"\s*\.le\.\s*" => " <= ",
  r"\s*\.ge\.\s*" => " >= ",
  r"\s*\.gt\.\s*" => " > ",
  r"\s*\.lt\.\s*" => " < ",
  # Remove (expression) brackets after if
  # r"if \((.*)\)(\s*\n)" => s"if \1\2",
  # Add end after single line if with an = assignment
  r"if\s*(.*?) = (.*?)(\n)" => s"if \1 = \2 end\3",
  # Format floats as "5.0" not "5."
  r"(\W\d+)\.(\D)" => s"\1.0\2",
  # Tab to 4 spaces
  r"\t" => "    ",
  # Relace suberror with error and mark for fixup
  r"(\W)suberror\((.*?),.*?\)" => s"\1 error(\2)",
  # Mark #FIXME the various things this script can't handle
  r"(write|goto|while\s)" => s"#FIXME \1",
)

# Patterns to remove
removal = [
  # Trailing whitespace
  r"\h*$"m,
  # Variable declarations
  r"\n\s*real\s.*",
  r"\n\s*real, external\s.*",
  r"\n\s*integer\s.*",
  r"\n\s*implicit none",
  r"\n\s*logical\s.*",
  # Import statements
  r"\n\s*use\s.*",
]

# Load the file from the first command line argument
filename = string(ARGS[1])
global code = read(filename, String)

# Process replacements and removals.
for (f, r) in replacements
  global code = replace(code, f, r)
end
for r in removal
  global code = replace(code, r, "")
end
println(code)


# Write the output to a .jl file with the same filename stem.
stem = split(filename, ".")[1]
outfile = stem * ".jl"
write(outfile, code)


#=
Copyright (c) 2022 Rafael Schouten
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#