# Graphics — plotting in Quanty

How Quanty draws: the `Graphics` object, the functions that make one, colour
maps, band structures with a self energy, output formats and traps. For what a
function does from the user's side, read the Graphics pages of the language
reference on <https://www.quanty.org>
(`documentation:language_reference:objects:graphics`). The wiki and the
option names write "color"; this note writes "colour" in prose and keeps option
names exactly as in the code, misspellings included (`Transparancy`,
`WatterfallShift`).

The colour convention, `0xAARRGGBB` with `AA` the transparency, is §5 of
`CONVENTIONS.md` and is not repeated here. Where `Graphics` and `BitMap` are
registered, and the names behind their `__index`, is in `LUA_API_OVERVIEW.md`.


## 1. Where the code is

| File | What |
|---|---|
| `Plot/StructDef.cpp` | `GraphicsType` and its elements, the colour helpers, the option structs `FramePlotOptions` and `LinePlotOptions`, the font list `AvailableFontTypes` |
| `Plot/Plot.cpp` | init, copy and free; `AddLine`, `AddCurve`, `AddArc`, `AddText`, `AddBitMap`, `AddLayer`; `RotateGraphics` |
| `Plot/ScientificPlots.cpp` | frames, ticks, labels, curves from spectra, line thinning, colour maps, band bitmaps |
| `Plot/GraphicsToSVG.cpp`, `Plot/GraphicsToPS.cpp`, `Plot/BitMap.cpp` | the writers: SVG, PostScript, BMP |
| `Plot/LuaBinding.cpp` | every Lua function: `Graphics_func[]`, `BitMap_func[]`, `LuaPlotGraphics`, `LuaPlotBands` |
| `LuaInterface/TightBinding.cpp` | `LuaCalculateBands` and its self-energy routes (§7) |

`LuaInterface/Plot.cpp` and `LuaInterface/PostScript.cpp` are not compiled (their
`#include` in `main.cpp` is commented out), so there is no `PS` table in the
binary; do not cite them as a precedent. The option parser for `Frame` and `LineStyle` lives with all other
options, in `LuaCheckOptions` (`LuaInterface/lua_extensions/TableToArray.cpp`).
Option names are compared without regard to case (`strcicmp`), so `plotRe` and
`PlotRe` are the same option.


## 2. The object

A `GraphicsType` is a canvas with lists of elements and a list of layers.

- **Canvas.** `CanvasSize` and `CanvasOrigin` in mm, x to the right, y up. The
  default (`InitGraphics`) is A4 landscape at the origin.
- **Positions** of every element are fractions of the canvas, 0 to 1. **Widths,
  text sizes and radii** are fractions of `sqrt(width * height)` of the canvas,
  so they keep their look when the canvas changes shape. **Angles** are radians,
  counter-clockwise from the x axis.
- **Elements:** `Line` (a polyline, optionally closed and filled), `Curve` (cubic
  Bezier: each point has a left and a right control point, `VL` and `VR`), `Arc`
  (centre, radius, two angles, optionally closed through the centre and filled),
  `Text` (string, point, angle, alignment `l`, `c` or `r`, size, font) and
  `BitMap` (pixels in the rectangle `LB`–`RT`, row 0 at the top, a `Pixelated`
  flag and an `Angle`).
- **Layers.** `Layer[0]` is the graphics itself; `Layer[1…NLayer-1]` are other
  graphics with their own canvas in the same absolute mm coordinates of the page.
  A layer is a copy: adding a graphics as a layer copies it (`AddLayer`).

**Draw order**, the same in SVG and PostScript (`AddLayersToSVG`,
`AddLayersToPS`): within one graphics all bitmaps, then lines, curves, arcs and
texts; then its layers in order, each recursively. So in `G1 + G2` a bitmap of
`G2` covers the lines of `G1`.

Moving or resizing a graphics (`SizeX`, `SizeY`, `OriginX`, `OriginY`, in
`LuaGraphicsNewIndex`) moves and scales its layers with it. `NLayer` can be read
as a property, besides the counts listed in `LUA_API_OVERVIEW.md`.


## 3. Building a graphics by hand

- **`Graphics.New()`** or `Graphics.New({{"Canvas.Size.X",w},{"Canvas.Size.Y",h},
  {"Canvas.Origin.X",x},{"Canvas.Origin.Y",y}})`. The option `NLayer` was removed
  in `ab671f03`; it was never used. `Graphics.New(T)` with a table in the form of
  `Graphics.ToTable` rebuilds the object; `io.put` relies on the pair.
- **`Graphics.ToTable(G)`** gives every element as plain Lua tables; colours as
  numbers, bitmap pixels as one hex string of eight digits per pixel.
- **`G.Add(...)` / `Graphics.Add(G, ...)`** (`LuaAddToGraphics`). The method
  changes `G`; the static form returns a changed copy and leaves `G` alone, the
  convention of §4 of `LUA_API_OVERVIEW.md`. Several objects can follow each
  other in one call, and an option list comes last:
  - `"text", x, y` — a text;
  - `{{x1,y1},{x2,y2},...}` — a line;
  - `{{vlx,vly,px,py,vrx,vry},...}` — a Bezier curve;
  - `{"Arc", {x,y}, R}` or `{"Arc", {x,y}, R, Theta1, Theta2}` — a circle or an
    arc; both angles are taken modulo 2 pi, so `{-Pi/2, Pi/2}` is a half circle;
  - `{"BitMap", LB, RT, Pixels}` or `{"BitMap", LB, RT, Pixels, TransparentColor}`
    — a table of rows (the first at the top) of pixels `grey`, `{grey}`,
    `{r,g,b}` or `{r,g,b,t}`, 0 to 255;
  - `{"BitMap", LB, RT, f, M1, M2, ...}` — pixel `[i][j]` gets the colour
    `f(Re M1[i][j], Im M1[i][j], Re M2[i][j], ...)`, which returns `r, g, b`
    from 0 to 1 (clipped). These are drawn with hard pixels by default;
  - a `Graphics` object (added as a layer) or a `BitMap` object.

  Options: `Color`, `Fill Color`, `Close`, `Line Width` (default 0.001), `Size`
  (text, default 0.05), `Alignment`, `Angle`, `Font` (default
  `TimesNewRomanPSMT`), `Transparancy`, `Fill Transparancy`, `alpha`,
  `Fill alpha`, `Pixelated`, and `Layer`: `n` adds to layer `n` (0 the graphics
  itself, -1 the last layer), `{n1, n2}` to layer `n2` of layer `n1`.
- **`G.Rotate(Angle [, opts])` / `Graphics.Rotate(G, Angle [, opts])`** turn a
  graphics, or one layer with `{"Layer", n}`, about `{"Center",{x,y}}` (default
  the middle). Everything keeps its size in mm; each canvas becomes the bounding
  box of its rotated self. A bitmap turned by a multiple of 90° has its pixels
  moved; by any other angle it keeps them and gets `BitMap.Angle`.
- **`G1 + G2`** is a copy of `G1` with `G2` as its last layer. `#G` is a stub
  that returns nothing.
- **`BitMap.New(nx, ny [, opts])`**, **`BitMap.Pixel(B, column, row, colour
  [, alpha])`** (1-based, row 1 at the top) and **`BitMap.Save(B, file)`** (a BMP
  file). A `BitMap` has no methods; `BitMap.FromMatrix`, `BitMap.Overlay`, `+`
  and `#` on it are empty stubs.


## 4. `Graphics.Plot`

`LuaPlotGraphics` takes one thing to plot and an option list, and returns a new
graphics with frame, ticks, labels and curves. What it accepts:

- a `Spectra` object, or a table of single-spectrum `Spectra` objects;
- a list of points `{{x1,y1},{x2,y2},...}` (y real or complex), or a table of
  such lists; the points may lie on any mesh;
- an `InterpolatingFunction`, or a table of them. It is evaluated on
  `PlotPoints + 1` points and on its own points inside the range. One with a
  continuous derivative is drawn as Bezier curves (`a6fa487a`);
- a response function, userdata or table: `Graphics.Plot(G [, GammaL [, GammaG
  [, NE]]], opts)`. Without `GammaL` or `GammaG` every pole is a spike on the
  grid, which is rarely wanted;
- a table of response functions with named broadenings,
  `Graphics.Plot({G1, G2, GammaL = 0.1, GammaG = 0}, opts)`. Each element gives
  one curve: for a block response function only its first spectrum, the (1,1)
  element.

A `Matrix` gives an error that names a function `MatrixPlot`, which does not
exist; a `GaussianRadialFunction` cannot be plotted yet.

**Options** (all in `LuaGraphicsOptions`): `Frame` (alias `Axis`), `LineStyle`,
`LineStyleCycle`, `plotRe`, `plotIm`, `Negate`, `Scale`, `PlotPoints` (600),
`Nrow`, `Ncolumn`, `ScaleOffDiagonal`, `ShiftOffDiagonal`, `WaterfallOffset`,
`WaterfallShift` (old spelling `WatterfallShift` still accepted),
`UseResonantEnergy`, `SpectraStart`, `SpectraStop`, `SpectraStep`, and the colour
map options of §5. `Negate` draws `-Im G`, the positive density of states.

**`Frame`** options: `Xmin`, `Xmax`, `Ymin`, `Ymax`; margins `XLMargin`,
`XRMargin`, `YLMargin` (bottom), `YRMargin` (top); `LineWidth`, `Color`,
`FontSize`, `Font`, `XLabel`, `YLabel`; `FirstXTick`, `dXTick`, `NMinorXTick`,
`FirstYTick`, `dYTick`, `NMinorYTick` (ticks 1 apart with 9 minor ticks by
default); `XFormat`, `YFormat` (default `"%.0f"`); `XAxis`, `YAxis`; `XTicks`,
`YTicks`. The x range defaults to that of the first spectrum, or -1 to 1 for a
response function. The y range is -1 to 1, except that a waterfall (a shift, an
offset or resonant energies) and a colour map choose it from the data.

**Line styles.** `LineStyle` is a list of `{{"LineWidth",w},{"LineColor",c},
{"PointSize",s},{"PointColor",c}}`, numbered from 1. The defaults are black, red
and blue, width 0.002, no points. Which style a curve gets:

- one spectrum with only Re or only Im: style 1; with both: style 2 for Re and 3
  for Im;
- several spectra with only Re or only Im: spectrum `i` gets style `i`, cycling
  after `LineStyleCycle` styles. With `plotRe = false` the list therefore starts
  at the first Im curve; no placeholder `{}` is needed;
- several spectra with both: Re of spectrum `i` gets style `2i-1`, Im style `2i`.

**Block response functions.** A block response function drawn directly gives a
`BlockSize × BlockSize` grid of panels, each a layer
(`PlotSpectraTypeArrayToGraphicsArray`); `ScaleOffDiagonal` and
`ShiftOffDiagonal` rescale the off-diagonal panels. To draw chosen elements in one
panel, see §6.

**Lines are thinned.** `AddDataPointsToGraphics` drops points that lie within
0.01 mm of the line through the kept ones (Douglas–Peucker, `435896ed`); the
SVG of a dense curve shrinks tenfold with the same picture. Dots (`PointSize`)
are still drawn at every point, which is how the tests count points.

**Lines are cut at the frame.** A part of a curve outside the plot range is left
out, and each part inside ends where the curve crosses the edge
(`FindBoundaryPoint`, or `FindBoundaryPointOfCubic` for Bezier pieces). The edge
is the one the segment moves towards; until this was fixed, a line whose point
lay on the left or bottom edge and that left the frame through the top shrank to
that point. Dots outside the range are not drawn.


## 5. Colour maps and bitmaps

`{"ColorMap", true}` draws a list of spectra as one bitmap instead of curves
(`PlotSpectraTypeArrayToColorMap`), for example a RIXS map. Spectrum `i` lies at
its waterfall height: its resonant energy with `UseResonantEnergy`, else
`WaterfallOffset + (i-1) WaterfallShift`, or without a shift `WaterfallOffset`
plus its number in the input. The heights must increase or decrease strictly; a map needs at least
two spectra and cannot be combined with `Nrow`/`Ncolumn` or with point lists.

- **The built-in scale** goes from white to `Color` as `Saturation × Im` goes
  from 0 to 1 (`Re` when only `plotRe` is set). `Scale` and `Negate` act before
  the colour, so `-Im` needs `Negate`.
- **`ColorFunction`** replaces it: `f(Re, Im)` returns `r, g, b` from 0 to 1,
  clipped, as the character functions of `TightBinding.CalculateBands` do (not 0
  to 255). It is stored as a registry reference; an error inside it is raised
  again after the pixels are freed, with its own message.
- `NPixelX`, `NPixelY` (600 each) set the resolution. `Interpolate` (default
  true) interpolates linearly in energy and between spectra; `false` takes the
  nearest value and gives `image-rendering: pixelated` in SVG.
- Pixels without data are fully transparent.

**There is no colour bar** yet; it was the agreed next step after the colour map.

**Bitmaps are written as BMP with a `BITMAPV5HEADER` and `BI_BITFIELDS`**
(`Plot/BitMap.cpp`), base64 inside the SVG. The old 40-byte header made Firefox
and WebKit ignore the alpha channel, and a V4 header did not decode in WebKit;
V5 works in both.


## 6. Plotting single elements of a block response function

To plot chosen elements of a block response function in one graph, go through
`Spectra`:

    S = ResponseFunction.ToSpectra(G, {{"Emin",-5},{"Emax",5},{"GammaL",0.1}})
    pl = Graphics.Plot(S.Element({1,4}), {{"plotRe",false},{"Negate",true},
                       {"LineStyle",{style1, style2}}})

`ToSpectra` flattens the block into `BlockSize^2` spectra; for a 2 × 2 block,
spectra 1 and 4 are the diagonal elements. Give one `LineStyle` entry per
spectrum. `Graphics.Plot` with `SpectraStart`, `SpectraStop` and `SpectraStep`
on a block response function selects in the same flattened numbering.

**Only principal sub-blocks are response functions.** An off-diagonal element
such as `G[1][2]` is not one: its imaginary part is not negative semi-definite.
A sub-block with the same index set for rows and columns is. So never select an
off-diagonal element as a `ResponseFunction`; plot it as a spectrum.
`ResponseFunction.Elements` is an empty stub; if it is ever written, it may only
return principal sub-blocks.


## 7. Band structures

`Graphics.PlotDOS`, `Graphics.PlotBandsAndDOS` and `Graphics.PlotFermi` are
stubs that return nothing, and have no wiki page.

**`Graphics.PlotBands(kpath, E [, c1 [, c2, c3]], opts)`** (`LuaPlotBands`)
draws bands as a bitmap inside a frame with the k-path labels. The k-path is
`{{Label, k}, N1, {Label, k}, N2, ...}`. **A label is one character**: "Gamma"
is drawn "G".

- With matrices from `TightBinding.CalculateBands` (energies, characters): no
  character draws the bands in `Color`; one character mixes `Color` and
  `SecondColor`; three are used as red, green and blue. Two or four characters
  (a band thickness) are not implemented: the thickness is ignored.
- With tables of one response function per k-point (from `CalculateBands` with a
  self energy) the spectral functions are drawn, broadened by `GammaL` and
  `GammaG`. The colour comes from `c = A_c / A_tot`, the opacity is
  `Saturation × A_tot` over the peak of a unit pole, blended with white
  (`LuaBandSpectraToBitMap`, `ConvertBandSpectraToBitMap`). Matrices with
  `GammaL` or `GammaG` are drawn the same way, as poles of weight 1. A thickness
  is an error here.
- Options: `NPixelX`, `NPixelY` (1000 each), `Color` (dark blue `0x000080`),
  `Transparancy`, `SecondColor` (red), `SecondTransparancy`, `LineWidth`,
  `EmptyPathWidth`, `Frame`, `Interpolate`, `GammaL`, `GammaG`, `Saturation`.
  `Transparancy` multiplies with the colour's own `AA` as opacities
  (`16e136ee`).
- **Set `NPixelX` and `NPixelY`.** The default 1000 × 1000 bitmap makes an SVG
  of about 5 MB. These options were parsed into the wrong type until `550902fb`.

**`TightBinding.CalculateBands(HTB, kpath, Sigma [, f1 [, f2, f3]], opts)`**
with a self energy returns `Gtot, G1 [, G2, G3]`: tables of one scalar list of
poles per k-point. `Sigma` at position 3 may be any form `CalculateG` accepts.
**A Lua table at position 3 always means `Sigma`**, so a character given as a
matrix must be a `Matrix` object, not a table. A character function must be a
quadratic form; it is turned into a matrix by polarisation and checked
(`LuaCharacterToMatrix`).

- **The default route** (`LuaCalculateBandsSelfEnergy`) diagonalises `H(k)` plus
  the embedded self energy at every k. Exact, but it costs `Nk × N^3` with
  `N = n (1 + number of poles of Sigma)`. Reduce the poles first with
  `ResponseFunction.ReducePoles` (`Method` `"Number"` or `"Weight"`).
- **`{"Method","Mesh"}`** (`LuaCalculateBandsSelfEnergyMesh`,
  `TightBindingBandsPlusSelfEnergyOnMesh`) puts the poles on an energy grid
  (`Emin`/`Emax`/`NE` or `EnergyGrid`, and `Lambda`, default 2). Its cost does
  not grow with the number of poles. The grid spacing is the resolution and must
  be well below the broadening used to draw.

**The gauge.** `H(k)` and the self energy are added in the gauge of
`CalculateG`, constant phase within the unit cell (`ConstantPhaseInUnitCell`),
in which a k-independent self energy couples orbitals of the same cell `R`. That
is also the gauge in which `CalculateG` reads an inter-site self energy. The
eigenvectors are then rotated, `u_i -> e^{I k.r_i} u_i`, to the gauge of
`CalculateBands` without a self energy, so a character means the same with and
without one. The mesh route has no eigenvectors and rotates the character
instead, `D^dag P D` with `D = diag(e^{I k.r_i})`, which gives the same
`Tr[P A]`. A local self energy does not depend on the gauge.


## 8. Output

- **`Graphics.ToSVG(G [, {{"RelativeSize",true}}])`** returns the SVG as a string.
  `RelativeSize` sets width and height to 100 %, for embedding in a web page.
  Texts are escaped (`&`, `<`, `>`) since `c6fb618a`. Use SVG.
- **`Graphics.ToPS(G)`** returns PostScript. **PostScript does not handle
  transparency correctly**: partial transparency is ignored, fully transparent
  fills and texts are left out, but a fully transparent stroke is drawn opaque,
  and `colorimage` drops the alpha channel of bitmaps. This is accepted; the wiki
  says only that, and recommends SVG.
- **`Graphics.SVGHTMLViewer()`** returns two HTML pages as strings: a viewer for
  one SVG file and a player for a sequence of them.
- Write files with `io.open`. Make a folder with **`io.mkdir(path)`**, never
  `os.execute("mkdir ...")`, which differs between Windows and Unix. It makes the
  parents too, returns `true` if it made the folder and `false` if it existed.

**Checking a plot.** `rsvg-convert` renders embedded bitmaps blank; headless
Firefox draws everything. To check a band bitmap by numbers, decode its base64
BMP (32-bit BGRA): the weight is in the luminance, the alpha is uniform.


## 9. Tight-binding traps that show up in band plots

- **No automatic Hermitian conjugate.** `CreateHamiltonianFromTightBindingDef`
  adds each entry of `HTB.Hopping` once. Every bond must be listed for `+T` and
  for `-T`, or `H(k)` is not Hermitian, with no warning. Do not copy the hopping
  lists of `TestScripts/TightBindingDef/DOSandBands.Quanty`, which give `+T` only.
- **Units.** The default is `Absolute`: the third entry of a hopping is a
  Cartesian vector from atom 1 to atom 2 in its cell, and Quanty checks that it
  is a whole number of cells plus the offset between the atoms (to 1e-7 cells),
  so a typo is caught. With `"Relative"` in `HTB.Units` it is a triple of
  integers. The `Relative` bugs in `H(k)`, the connecting vector and the band
  gauge were fixed in `721f5db9`, supercells and clusters in `de2d12ce`
  (`AtomPositionCartesian`; test `RelativeUnits.Quanty`). What remains open is in
  `TODO.md`; `Absolute` stays the safer choice.
- **`H(k)` uses `TransVecs`**, never `ConnectingVector`, which several C
  builders leave at 0.
- **A cheap check on a hopping list:** `TightBinding.CalculateMu` (it exists only
  in that table) returns `mu`, `Emin`, `Emax` and the gap.


## 10. Rules for scripts and tutorials

- **One scale per figure.** Curves in one figure share their scale. If one is
  multiplied for readability, write the factor in the plot — in its label,
  legend or a note under the figure ("FY × 4") — not only in the script. A reader
  compares heights; a hidden factor misleads. Experimental data, on an arbitrary
  scale, may be scaled to the theory without a label.
- **Call `Graphics.Plot` directly**, with a `Frame`, labels, and `Graphics.ToSVG`.
  Do not wrap it in a new plotting function with its own automatic ranges, tick
  arithmetic or scaling: the scripts are there to show the Quanty function. If
  `Graphics.Plot` falls short, that is a C issue to report and fix, not to work
  around in Lua.
- **Broaden what you plot**: `Graphics.Plot({G, GammaL = ...}, {{"Negate",true},
  {"plotRe",false}})` for a density of states, not a bare `Graphics.Plot(G)`.
- New colours are `0xAARRGGBB` (`CONVENTIONS.md` §5); an old script needs
  converting.


## 11. Tests

`TestScripts/Plot/` (`Add`, `BitMap`, `New`, `Plot`, `Rotate`, `TextFont`,
`ToTable`) is the `test-jobs-plot` block of `.gitlab-ci.yml`. `Plot.Quanty`
reads plots back with `Graphics.ToTable`. `Add.Quanty` rewrites tracked SVGs
under `TestScripts/Plot/OutputPlots/`; a changed file there is a signal. Bands
with a self energy: `TestScripts/TightBindingDef/CalculateBandsSelfEnergy.Quanty`
(every input form of `PlotBands`) and `CalculateBandsMesh.Quanty`. Do not compare
the two routes above about 200 poles of `Sigma`; the default route then takes
very long.
