# Third-party notices

The Steam-generating heat pump web app (`apps/steam/`) is licensed under the Apache License,
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

## Google Fonts

The page loads Barlow Semi Condensed, Source Sans 3 and IBM Plex Mono from the Google Fonts
CDN (`fonts.googleapis.com` and `fonts.gstatic.com`) at run time. The font files are not
included in this repository. These fonts are licensed under the
[SIL Open Font License, Version 1.1](https://openfontlicense.org/).
