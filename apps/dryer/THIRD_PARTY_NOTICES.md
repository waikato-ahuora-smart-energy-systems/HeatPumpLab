# Third-party notices

The Advanced dryer heat pumps web app (`apps/dryer/`) is licensed under the Apache License,
Version 2.0 (see the [LICENSE](../../LICENSE) at the repository root). It bundles or uses
the third-party components below, which keep their own licences.

## CoolProp

`index.html` contains the Emscripten JavaScript glue for CoolProp, and `coolprop.wasm` is
CoolProp compiled to WebAssembly. Source: <https://github.com/CoolProp/CoolProp>.
The licence below is copied from the CoolProp repository's
[LICENSE](https://github.com/CoolProp/CoolProp/blob/master/LICENSE) file.

```
MIT License

Copyright (c) 2012-2018 Ian H. Bell and other CoolProp developers

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
```

## teqp (ammonia–water coefficients)

The ammonia–water properties in `index.html` (from `apps/src/dryer/ammonia_water.js`) embed the coefficients of the
IAPWS G4-01 (Tillner-Roth and Friend, 1998) mixture model and of its pure-fluid equations of state (Tillner-Roth et
al. 1993 for ammonia, IAPWS-95 for water) as hard-coded in teqp 0.23.2
(`include/teqp/models/ammonia_water.hpp`). Source: <https://github.com/usnistgov/teqp>. The notice below is copied
from teqp's [LICENSE.md](https://github.com/usnistgov/teqp/blob/main/LICENSE.md).

```
NIST Disclaimer of Copyright and Warranty

This software was developed by employees of the National Institute of
Standards and Technology (NIST), an agency of the Federal Government
and is being made available as a public service. Pursuant to title 17
United States Code Section 105, works of NIST employees are not
subject to copyright protection in the United States. This software
may be subject to foreign copyright. Permission in the United States
and in foreign countries, to the extent that NIST may hold copyright,
to use, copy, modify, create derivative works, and distribute this
software and its documentation without fee is hereby granted on a
non-exclusive basis, provided that this notice and disclaimer of
warranty appears in all copies.

THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND,
EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED
TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS,
ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE
DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE
SOFTWARE WILL BE ERROR FREE. IN NO EVENT SHALL NIST BE LIABLE FOR ANY
DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR
CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY WAY
CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY,
CONTRACT, TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY
PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS
SUSTAINED FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE
SOFTWARE OR SERVICES PROVIDED HEREUNDER.
```

## Google Fonts

The page loads Barlow Semi Condensed, Source Sans 3 and IBM Plex Mono from the Google Fonts
CDN (`fonts.googleapis.com` and `fonts.gstatic.com`) at run time. The font files are not
included in this repository. These fonts are licensed under the
[SIL Open Font License, Version 1.1](https://openfontlicense.org/).
